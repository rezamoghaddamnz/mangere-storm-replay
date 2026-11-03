"""Self-checks for lisim: still water stays still, volume closes, and a rained-on plane reaches
Manning normal depth."""
import numpy as np

from lisim import Model


def test_lake_at_rest():
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 0.8, (40, 50))
    m = Model(z, 2.0, 0.03)
    m.h = np.maximum(1.0 - z, 0.0)
    h0 = m.h.copy()
    for _ in range(500):
        m.step(0.0)
    assert np.abs(m.qx).max() < 1e-12 and np.abs(m.qy).max() < 1e-12
    assert np.allclose(m.h, h0)


def test_volume_closes_with_rain_sinks_and_sea():
    rng = np.random.default_rng(2)
    z = rng.uniform(0, 2, (60, 60)) + np.linspace(3, 0, 60)[None, :]
    sea = np.zeros_like(z, bool)
    sea[:, -3:] = True
    m = Model(z, 4.0, 0.04, sea_mask=sea, sea_level=0.5)
    m.set_inlets([10, 30, 50], [20, 25, 40], [1e-3, 1e-3, 2e-3], [0, 0, 1], [0.02, 0.05])
    while m.t < 1800:
        m.step(50e-3 / 3600)
    while m.t < 3600:
        m.step(0.0)
    err = m.volume() - (m.rain_in - m.captured.sum() - m.sea_exchange)
    assert abs(err) / m.rain_in < 1e-9, err / m.rain_in
    assert m.captured[0] <= 0.02 * m.t + 1e-9 and m.captured[1] <= 0.05 * m.t + 1e-9


def test_volume_closes_with_moving_tide():
    rng = np.random.default_rng(3)
    z = rng.uniform(0, 1, (50, 50)) + np.linspace(2, -1, 50)[None, :]
    sea = z < -0.5
    m = Model(z, 4.0, 0.04, sea_mask=sea, sea_level=0.0)
    while m.t < 6 * 3600:
        m.sea_level = 1.2 * np.sin(2 * np.pi * m.t / (12.4 * 3600))
        m.step(20e-3 / 3600 if m.t < 3600 else 0.0)
    err = m.volume() - (m.rain_in - m.sea_exchange)
    assert abs(err) / m.rain_in < 1e-9, err / m.rain_in


def test_steep_crater_fills_without_spikes():
    """A 40 m deep cone with 35-degree sides: rain settles into a level pond of the right depth."""
    n, dx = 61, 2.0
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - n // 2, yy - n // 2) * dx
    z = np.minimum(r * np.tan(np.radians(35)), 40.0)
    m = Model(z, dx, 0.035)
    peak = 0.0
    while m.t < 1800:
        m.step(100e-3 / 3600)
        peak = max(peak, float(m.h.max()))
    while m.t < 3 * 3600:
        m.step(0.0)
        peak = max(peak, float(m.h.max()))
    vol = m.volume()
    level = (3 * vol / (np.pi / np.tan(np.radians(35)) ** 2)) ** (1 / 3)   # depth of a level pond in the cone
    assert abs(float(m.h.max()) / level - 1) < 0.1, (m.h.max(), level)
    assert peak < 1.2 * level, (peak, level)


def test_plane_reaches_normal_depth():
    nx, dx, S, n, r = 200, 2.0, 0.01, 0.03, 100e-3 / 3600
    z = np.tile(S * dx * np.arange(nx)[::-1], (3, 1))
    sea = np.zeros_like(z, bool)
    sea[:, -1] = True
    m = Model(z, dx, n, sea_mask=sea, sea_level=z[0, -1])
    while m.t < 4 * 3600:
        m.step(r)
    i = nx // 2
    x = (i + 0.5) * dx
    h_exact = (r * x * n / S ** 0.5) ** 0.6
    assert abs(m.h[1, i] / h_exact - 1) < 0.03, (m.h[1, i], h_exact)


if __name__ == "__main__":
    for f in (test_lake_at_rest, test_volume_closes_with_rain_sinks_and_sea, test_volume_closes_with_moving_tide,
              test_steep_crater_fills_without_spikes, test_plane_reaches_normal_depth):
        f()
        print("ok", f.__name__)
