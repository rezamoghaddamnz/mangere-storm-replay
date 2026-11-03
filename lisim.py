"""Local-inertial 2D overland flow (Bates et al. 2010; de Almeida et al. 2012), written from the
published scheme, with capacity-limited pipe inlets as sinks.

Grid: cell-centred depths h[ny, nx] over bed z; unit discharges qx on the ny x (nx+1) vertical
faces and qy on the (ny+1) x nx horizontal faces, in m2/s. Outer faces are closed walls; the sea
enters as cells held at a fixed water level.

Pipes are not routed hydraulically. Each inlet cell may remove water at its own inlet rate, and
all inlets draining to one outfall share that outfall pipe's full-bore capacity. That is a
capacity bound on the pipe system, not a pipe model, and every output that uses it says so.
"""
import numpy as np
from numba import njit, prange

G = 9.81
FR_MAX = 1.0
THETA = 0.7          # de Almeida et al. 2012: weights in the neighbouring fluxes, damps deep-water oscillation


@njit(parallel=True, cache=True)
def _flux_x(h, z, qx, q0, n, dt, dx, hmin):
    ny, nx = h.shape
    for j in prange(ny):
        for i in range(1, nx):
            eta_l, eta_r = h[j, i - 1] + z[j, i - 1], h[j, i] + z[j, i]
            hf = max(eta_l, eta_r) - max(z[j, i - 1], z[j, i])
            if hf <= hmin:
                qx[j, i] = 0.0
                continue
            nf = 0.5 * (n[j, i - 1] + n[j, i])
            q = q0[j, i]
            qa = THETA * q + 0.5 * (1.0 - THETA) * (q0[j, i - 1] + q0[j, i + 1]) if 1 < i < nx - 1 else q
            qn = (qa - G * hf * dt * (eta_r - eta_l) / dx) / (1.0 + G * dt * nf * nf * abs(q) / hf ** (7.0 / 3.0))
            qmax = FR_MAX * hf * (G * hf) ** 0.5   # Froude limiter: the scheme is unstable on steep slopes above Fr 1
            qx[j, i] = min(max(qn, -qmax), qmax)


@njit(parallel=True, cache=True)
def _flux_y(h, z, qy, q0, n, dt, dx, hmin):
    ny, nx = h.shape
    for j in prange(1, ny):
        for i in range(nx):
            eta_u, eta_d = h[j - 1, i] + z[j - 1, i], h[j, i] + z[j, i]
            hf = max(eta_u, eta_d) - max(z[j - 1, i], z[j, i])
            if hf <= hmin:
                qy[j, i] = 0.0
                continue
            nf = 0.5 * (n[j - 1, i] + n[j, i])
            q = q0[j, i]
            qa = THETA * q + 0.5 * (1.0 - THETA) * (q0[j - 1, i] + q0[j + 1, i]) if 1 < j < ny - 1 else q
            qn = (qa - G * hf * dt * (eta_d - eta_u) / dx) / (1.0 + G * dt * nf * nf * abs(q) / hf ** (7.0 / 3.0))
            qmax = FR_MAX * hf * (G * hf) ** 0.5   # Froude limiter: the scheme is unstable on steep slopes above Fr 1
            qy[j, i] = min(max(qn, -qmax), qmax)


@njit(parallel=True, cache=True)
def _limit_outflow(h, qx, qy, dt, dx):
    """Scale a cell's outgoing fluxes so it cannot give away more water than it holds."""
    ny, nx = h.shape
    for j in prange(ny):
        for i in range(nx):
            out = (max(qx[j, i + 1], 0.0) + max(-qx[j, i], 0.0)
                   + max(qy[j + 1, i], 0.0) + max(-qy[j, i], 0.0)) * dt / dx
            if out > h[j, i] and out > 0.0:
                f = h[j, i] / out
                if qx[j, i + 1] > 0.0:
                    qx[j, i + 1] *= f
                if qx[j, i] < 0.0:
                    qx[j, i] *= f
                if qy[j + 1, i] > 0.0:
                    qy[j + 1, i] *= f
                if qy[j, i] < 0.0:
                    qy[j, i] *= f


@njit(parallel=True, cache=True)
def _continuity(h, qx, qy, dt, dx, rain_depth):
    ny, nx = h.shape
    for j in prange(ny):
        for i in range(nx):
            dh = (qx[j, i] - qx[j, i + 1] + qy[j, i] - qy[j + 1, i]) * dt / dx
            h[j, i] = max(h[j, i] + dh + rain_depth, 0.0)


@njit(cache=True)
def _sinks(h, inlet_j, inlet_i, inlet_rate, inlet_group, group_cap, dt, cell_area, captured):
    """Inlets take min(water present, inlet rate); a group never exceeds its outfall capacity."""
    ng = group_cap.shape[0]
    demand = np.zeros(ng)
    want = np.empty(inlet_j.shape[0])
    for k in range(inlet_j.shape[0]):
        w = min(h[inlet_j[k], inlet_i[k]], inlet_rate[k] * dt)
        want[k] = w
        demand[inlet_group[k]] += w * cell_area
    scale = np.ones(ng)
    for g in range(ng):
        lim = group_cap[g] * dt
        if demand[g] > lim and demand[g] > 0.0:
            scale[g] = lim / demand[g]
    for k in range(inlet_j.shape[0]):
        take = want[k] * scale[inlet_group[k]]
        h[inlet_j[k], inlet_i[k]] -= take
        captured[inlet_group[k]] += take * cell_area


class Model:
    def __init__(self, z, dx, n, sea_mask=None, sea_level=0.0, hmin=1e-3, alpha=0.35):  # 0.7 went unstable in the Māngere Mountain crater; 0.35 and 0.2 agree
        self.z = np.ascontiguousarray(z, dtype=np.float64)
        self.dx, self.hmin, self.alpha = float(dx), hmin, alpha
        self.n = np.broadcast_to(np.asarray(n, dtype=np.float64), self.z.shape).copy()
        ny, nx = self.z.shape
        self.h = np.zeros((ny, nx))
        self.qx = np.zeros((ny, nx + 1))
        self.qy = np.zeros((ny + 1, nx))
        self.sea = np.zeros((ny, nx), dtype=np.bool_) if sea_mask is None else sea_mask.astype(np.bool_)
        self.sea_level = sea_level
        self.t = 0.0
        self.rain_in = 0.0            # m3 of rain added on land cells
        self.sea_exchange = 0.0       # m3 removed (+) or added (-) by holding the sea level
        self.inlets = None
        self._hold_sea()
        self.sea_exchange = 0.0       # filling the sea at start is not an exchange with land

    def set_inlets(self, j, i, rate, group, group_cap):
        """rate in m/s of depth per inlet cell; group_cap in m3/s per outfall."""
        self.inlets = (np.asarray(j, np.int64), np.asarray(i, np.int64), np.asarray(rate, np.float64),
                       np.asarray(group, np.int64), np.asarray(group_cap, np.float64))
        self.captured = np.zeros(len(group_cap))

    def _hold_sea(self):
        target = np.maximum(self.sea_level - self.z, 0.0)
        self.sea_exchange += float(((self.h - target) * self.sea).sum()) * self.dx ** 2
        self.h[self.sea] = target[self.sea]

    def dt(self):
        hmax = max(float(self.h.max()), 0.01)
        return self.alpha * self.dx / np.sqrt(G * hmax)

    def step(self, rain_rate, dt_max=5.0):
        """Advance one stable step; rain_rate in m/s on every non-sea cell."""
        dt = min(self.dt(), dt_max)
        self.h[self.sea] = np.maximum(self.sea_level - self.z, 0.0)[self.sea]   # tide moves: not an exchange with land
        _flux_x(self.h, self.z, self.qx, self.qx.copy(), self.n, dt, self.dx, self.hmin)
        _flux_y(self.h, self.z, self.qy, self.qy.copy(), self.n, dt, self.dx, self.hmin)
        _limit_outflow(self.h, self.qx, self.qy, dt, self.dx)
        _continuity(self.h, self.qx, self.qy, dt, self.dx, 0.0)
        if rain_rate > 0.0:
            land = ~self.sea
            self.h[land] += rain_rate * dt
            self.rain_in += rain_rate * dt * float(land.sum()) * self.dx ** 2
        if self.inlets is not None:
            j, i, rate, grp, cap = self.inlets
            _sinks(self.h, j, i, rate, grp, cap, dt, self.dx ** 2, self.captured)
        self._hold_sea()
        self.t += dt
        return dt

    def volume(self):
        return float((self.h * ~self.sea).sum()) * self.dx ** 2
