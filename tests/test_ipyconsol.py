import numpy as np
import pytest
import ucla_geotech_tools.ipyconsol as pcl

TMAX = 35532000.0
UNIFORM = dict(
    H=4.35,
    N=100,
    Cc=5.729,
    Cr=0.532,
    sigvref=50.0,
    esigvref=8.763,
    Gs=1.686,
    kref=1.0e-7,
    ekref=8.344,
    Ck=200.0,
    tref=3070.0,
    qo=10.0,
    ocrvoidratiotype=0,
    ocrvoidratio=1.5,
    ru=0.0,
    gammaw=9.81,
    pa=101.325,
    tol=1.0e-8,
)


def layered_inputs():
    H1 = 2.0
    H2 = 2.35
    N1 = 50
    N2 = 50
    depth1 = np.linspace(0, H1, num=N1, dtype=float, endpoint=False)
    depth2 = np.linspace(H1, H1 + H2, num=N2 + 1, dtype=float)
    depth = np.hstack((depth1, depth2))
    N = depth.size - 1
    return dict(
        depth=depth,
        Cc=np.full(N + 1, 5.729, dtype=float),
        Cr=np.full(N + 1, 0.532, dtype=float),
        sigvref=np.full(N + 1, 50, dtype=float),
        esigvref=np.full(N + 1, 8.763, dtype=float),
        Gs=np.full(N + 1, 1.868, dtype=float),
        kref=np.full(N + 1, 1.0e-7, dtype=float),
        ekref=np.full(N + 1, 8.344, dtype=float),
        Ck=np.full(N + 1, 200, dtype=float),
        Ca=np.full(N + 1, 0.325, dtype=float),
        tref=np.full(N + 1, 3070, dtype=float),
        dsigv=np.full(N + 1, 27, dtype=float),
        ocrvoidratiotype=np.hstack((np.full(N1, 0, dtype=int), np.full(N2 + 1, 1, dtype=int))),
        ocrvoidratio=np.hstack((np.full(N1, 1.5, dtype=float), np.full(N2 + 1, 10.0, dtype=float))),
        ru=np.full(N + 1, 0, dtype=float),
        qo=10,
        gammaw=9.81,
        pa=101.325,
        tol=1.0e-8,
    )


def ramp_load(Ntime=100):
    time = np.logspace(np.log10(TMAX) - 5, np.log10(TMAX), num=Ntime)
    return time, np.minimum(time / 100000, 1.0)


def check_output(output, Ntime):
    """Every time step converged and was stored, and the bottom of the layer does not move"""
    for key in ["z", "u", "sigv", "e"]:
        assert output[key].shape[1] == Ntime
        assert np.all(np.isfinite(output[key]))
    assert np.all(output["e"] > 0)
    assert np.all(output["sigv"] > 0)
    assert np.allclose(output["z"][-1], output["z"][-1, 0])


# test 1, uniform soil constant load
@pytest.mark.parametrize("drainagetype", [0, 1, 2])
def test_ipyconsol1(drainagetype):
    Ntime = 500
    output = pcl.compute(**UNIFORM, Ca=0.0, dsigv=400.0, Ntime=Ntime, tmax=TMAX, drainagetype=drainagetype)
    check_output(output, Ntime)
    # Fully consolidated at the end: effective stress increased by dsigv and void ratio on the normal consolidation line
    inputs = pcl.get_inputs(**UNIFORM, Ca=0.0, dsigv=400.0, Ntime=Ntime, tmax=TMAX)
    initial = pcl.get_initial(*inputs)
    assert np.allclose(output["u"][:, -1], 0.0, atol=1.0e-3)
    assert np.allclose(output["sigv"][:, -1], initial["sigvo"] + 400.0, atol=1.0e-3)
    e_ncl = UNIFORM["esigvref"] - UNIFORM["Cc"] * np.log10(output["sigv"][:, -1] / UNIFORM["sigvref"])
    assert np.allclose(output["e"][:, -1], e_ncl, atol=1.0e-6)


# test2 non-uniform soil, constant load
def test_ipyconsol2():
    Ntime = 100
    output = pcl.compute(**layered_inputs(), Ntime=Ntime, tmax=TMAX, drainagetype=0)
    check_output(output, Ntime)
    # Drained boundaries, and the surface settles monotonically
    assert np.all(output["u"][[0, -1]] == 0.0)
    assert np.all(np.diff(output["z"][0]) >= -1.0e-12)


# test 3 uniform soil, time-varying load
def test_ipyconsol3():
    time, loadfactor = ramp_load()
    output = pcl.compute(**UNIFORM, Ca=0.325, dsigv=27.0, time=time, loadfactor=loadfactor, drainagetype=0)
    check_output(output, len(time))
    assert np.all(output["u"][[0, -1]] == 0.0)
    assert np.all(np.diff(output["z"][0]) >= -1.0e-12)


# test 4 non-uniform soil, time-varying load
def test_ipyconsol4():
    time, loadfactor = ramp_load()
    output = pcl.compute(**layered_inputs(), time=time, loadfactor=loadfactor, drainagetype=1)
    check_output(output, len(time))
    assert np.all(output["u"][0] == 0.0)


def terzaghi_U(Tv):
    U = np.zeros_like(Tv)
    for m in range(200):
        M = np.pi / 2 * (2 * m + 1)
        U += 2 / M**2 * np.exp(-(M**2) * Tv)
    return 1 - U


@pytest.mark.parametrize("drainagetype", [0, 1, 2])
def test_terzaghi(drainagetype):
    # A small load increment and nearly constant permeability approach Terzaghi's linear theory.
    # Gs = 1 gives uniform initial effective stress.
    qo, dsigv, H, N, Ntime = 100.0, 1.0, 2.0, 100, 200
    Cc, eo, kref, gammaw = 0.5, 1.0, 1.0e-8, 9.81
    cv = kref * (1 + eo) / (Cc / np.log(10) / qo * gammaw)
    Hdr = H / 2 if drainagetype == 0 else H
    tmax = 2.0 * Hdr**2 / cv
    output = pcl.compute(
        Cc=Cc, Cr=Cc - 1.0e-9, sigvref=qo, esigvref=eo, Gs=1.0, kref=kref, ekref=eo, Ck=1.0e6, Ca=0.0, tref=1.0,
        qo=qo, dsigv=dsigv, ocrvoidratiotype=1, ocrvoidratio=eo, N=N, H=H, Ntime=Ntime, tmax=tmax,
        drainagetype=drainagetype, gammaw=gammaw,
    )
    time = np.logspace(np.log10(tmax) - 5, np.log10(tmax), Ntime)
    Uavg = 1 - np.trapezoid(output["u"], dx=H / N, axis=0) / (H * dsigv)
    assert np.max(np.abs(Uavg - terzaghi_U(cv * time / Hdr**2))) < 0.01


def test_scalar_functions_return_scalars():
    av = pcl.get_avtest(1.0, 100.0, 110.0, 1.0, 100.0, 0.5, 0.05)
    assert isinstance(av, float)
    assert np.isclose(av, 0.5 * np.log10(1.1) / 10.0)
    e = pcl.get_etest(1.0, 1.0, 1.0, 0.01, 1.0, 110.0, 100.0, 100.0, 0.5, av, 0.0, 10.0, 1.0)
    assert isinstance(e, float)
    assert np.isclose(pcl.get_ktest(1.0, 1.0, 1.0e-8, 0.5), 1.0e-8)
