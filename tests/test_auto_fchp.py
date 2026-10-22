import os

import numpy as np
import pytest
import ucla_geotech_tools.auto_fchp as af

dt = 0.02
acc = np.load(
    os.path.join(os.path.dirname(__file__), "2D.OBS03..BH1__20110107T211622Z__20110107T212922Z.npy")
)


def test_auto_fchp1():
    fchp = af.get_fchp(dt=dt, acc=acc)
    assert fchp == pytest.approx(0.0402, abs=0.001)


def test_auto_fchp2():
    fchp = af.get_fchp(
        dt=dt,
        acc=acc,
        target=0.02,
        tol=0.001,
        poly_order=6,
        maxiter=30,
        fchp_min=0.001,
        fchp_max=0.5,
        filter_order=5.0,
        tukey_alpha=0.05,
        apply_disp_ratio=1,
        disp_ratio_time=2,
        disp_ratio_target=0.02,
    )
    # Objective 2 can only raise fchp
    assert fchp == pytest.approx(0.0431, abs=0.001)
    assert fchp >= af.get_fchp(dt=dt, acc=acc)


def test_list_and_odd_length_input():
    assert af.get_fchp(dt=dt, acc=list(acc[:-1])) == pytest.approx(af.get_fchp(dt=dt, acc=acc[:-1]))


def test_invalid_keyword_returns_none():
    assert af.get_fchp(dt=dt, acc=acc, polyorder=6) is None


def test_solve_objective_bounds():
    # Satisfied everywhere: lowest value. Never satisfied: highest value. Otherwise the root.
    assert af.solve_objective(lambda f: -1.0, 0.01, 0.5, (), 1e-6, 30) == 0.01
    assert af.solve_objective(lambda f: 1.0, 0.01, 0.5, (), 1e-6, 30) == 0.5
    assert af.solve_objective(lambda f: 0.1 - f, 0.01, 0.5, (), 1e-9, 30) == pytest.approx(0.1)
    # Satisfied at the lower bound even if the residual is positive at the upper bound
    assert af.solve_objective(lambda f: f - 0.2, 0.01, 0.5, (), 1e-6, 30) == 0.01
