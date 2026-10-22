import numpy as np
import scipy as sp
import ucla_geotech_tools.response_spectrum as rs


def get_motions(Nmotion):
    N = 2**14
    dt = 0.005
    time = np.linspace(0, dt * N, N)
    f = np.linspace(0.5, 8.0, Nmotion)
    motions = np.sin(2 * np.pi * f[:, np.newaxis] * time)
    window = sp.signal.windows.tukey(N, 0.20)
    mean_motions = np.average(motions, weights=window, axis=1)
    motions = window * (motions - mean_motions[:, np.newaxis])
    return motions, dt


def test_response_spectrum1():
    motions, dt = get_motions(20)
    T = rs.get_ngawest2_T()
    Sa = rs.get_response_spectrum(motions=motions, D=0.05, T=T, dt=dt)
    assert Sa.shape == (20, len(T))
    # Interpolation only raises peaks that fall between time steps
    Sa1 = rs.get_response_spectrum(motions=motions, D=0.05, T=T, dt=dt, interpolation_factor=1)
    assert np.all(Sa >= Sa1 * (1 - 1e-3))


def test_pga():
    motions, dt = get_motions(3)
    Sa = rs.get_response_spectrum(motions=motions, T=[0.0], dt=dt)
    assert np.allclose(Sa[:, 0], np.max(np.abs(motions), axis=1))


def test_pseudo_spectral_acceleration():
    # Pseudo-spectral acceleration is omega_n^2 times peak relative displacement from a time-domain solution
    motions, dt = get_motions(2)
    T = np.asarray([0.1, 0.5, 2.0])
    D = 0.05
    Sa = rs.get_response_spectrum(motions=motions, D=D, T=T, dt=dt)
    time = dt * np.arange(motions.shape[1])
    for k, Tn in enumerate(T):
        wn = 2 * np.pi / Tn
        system = sp.signal.lti([-1.0], [1.0, 2 * D * wn, wn**2])
        for i, motion in enumerate(motions):
            _, u, _ = sp.signal.lsim(system, motion, time)
            assert np.isclose(Sa[i, k], wn**2 * np.max(np.abs(u)), rtol=0.01)


def test_response_spectrum2():
    motions, dt = get_motions(2)
    T = rs.get_ngawest2_T()
    Sa = rs.get_response_spectrum(motions=motions, D=0.05, T=T, dt=dt, RotD=50, N_rot_angle=100)
    assert Sa.shape == (len(T), 1)
    Sa = rs.get_response_spectrum(motions=motions, D=0.05, T=T, dt=dt, zeropad=1, RotD=[0, 50, 100])
    assert Sa.shape == (len(T), 3)
    assert np.all(np.diff(Sa, axis=1) >= 0)
    Sa_single = rs.get_response_spectrum(motions=motions, D=0.05, T=T, dt=dt)
    # RotD100 is at least the larger of the two as-recorded components
    assert np.all(Sa[:, 2] >= np.max(Sa_single, axis=0) * (1 - 1e-9))


def test_rotd_against_angles():
    rng = np.random.default_rng(0)
    data1 = rng.standard_normal(2000)
    data2 = 0.5 * rng.standard_normal(2000) + 0.3 * data1
    theta = np.linspace(0, np.pi, 3601)
    peaks = np.max(np.abs(np.outer(np.cos(theta), data1) + np.outer(np.sin(theta), data2)), axis=1)
    RotD = rs.rotd_percentiles(data1, data2, [0, 50, 100])
    assert np.allclose(RotD, np.percentile(peaks, [0, 50, 100]), rtol=1e-3)


def test_rotd_single_component():
    # With one component zero, the peak at angle theta is |cos(theta)| max|data1|
    data1 = np.sin(np.linspace(0, 20, 1000))
    RotD = rs.rotd_percentiles(data1, np.zeros_like(data1), [0, 50, 100])
    assert np.allclose(RotD, [0.0, np.max(np.abs(data1)) * np.cos(np.pi / 4), np.max(np.abs(data1))], atol=1e-9)


def test_interpolation_factor():
    assert rs.get_interpolation_factor(0.001) == 1
    assert rs.get_interpolation_factor(0.005) == 8
    assert rs.get_interpolation_factor(0.01) == 8
    assert rs.get_interpolation_factor(0.002) == 2


def test_sinc_interpolate():
    # Sinc interpolation of a band-limited signal reproduces the original samples and the signal between them
    dt = 0.01
    t = dt * np.arange(1024)
    data = np.sin(2 * np.pi * 3.0 * t) * sp.signal.windows.tukey(1024, 0.2)
    fine = rs.sinc_interpolate(data, 4)
    assert len(fine) == 4 * len(data)
    assert np.allclose(fine[::4], data, atol=1e-6)
    fine_t = dt / 4 * np.arange(len(fine))
    exact = np.sin(2 * np.pi * 3.0 * fine_t) * np.interp(fine_t, t, sp.signal.windows.tukey(1024, 0.2))
    assert np.allclose(fine[100:-100], exact[100:-100], atol=1e-3)
