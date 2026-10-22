import numpy as np
import pytest
import ucla_geotech_tools.random_field as rf


def test_random_field1():
    rng = np.random.default_rng(0)
    z = rng.random(1000)
    x = rng.random(1000)
    for model in (0, 1):
        field = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1, model=model, seed=1)
        assert field.shape == x.shape
        assert np.all(np.isfinite(field))


def test_random_field2():
    rng = np.random.default_rng(0)
    z = rng.random(1000)
    x = rng.random(1000) * 10
    theta_z = 0.5
    field2a = rf.get_random_field(x=x, z=z, theta_z=theta_z, mu=0, sigma=1, model=0, theta_x=10 * theta_z, seed=1)
    field2b = rf.get_random_field(
        x=x, z=z, theta_z=theta_z, mu=0, sigma=1, model=0, theta_x=10 * theta_z, Dhchol=0.1, seed=1
    )
    assert np.all(np.isfinite(field2a)) and np.all(np.isfinite(field2b))


@pytest.mark.parametrize("n", [10, 21, 41])
def test_partition_on_regular_grid(n):
    # Distances tie on regular grids; every node must still be assigned exactly once
    X, Z = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    x, z = X.ravel(), Z.ravel()
    Dh = 0.15 * 0.5
    chol, xc, zc, krig, xk, zk = rf.get_chol_krig_nodes(x, z, Dh)
    assert np.array_equal(np.sort(np.concatenate([chol, krig])), np.arange(n * n))
    # Every kriging node is within Dhchol of a cholesky node
    dist = np.hypot(xk[:, None] - xc, zk[:, None] - zc)
    assert np.all(dist.min(axis=1) <= Dh)
    field = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1, seed=0)
    assert np.all(np.isfinite(field))


@pytest.mark.parametrize("model", [0, 1])
def test_statistics(model):
    # Mean, variance, and correlation of many realizations match the target
    rng = np.random.default_rng(0)
    x, z = rng.random(150) * 2, rng.random(150) * 2
    theta, nreal = 0.5, 2000
    fields = np.array(
        [rf.get_random_field(x=x, z=z, theta_z=theta, mu=5.0, sigma=2.0, model=model, seed=k) for k in range(nreal)]
    )
    assert fields.mean() == pytest.approx(5.0, abs=0.05)
    cov = np.cov(fields.T) / 4.0
    assert np.diag(cov).mean() == pytest.approx(1.0, abs=0.03)
    target = rf.get_C(x, z, theta, model)
    standard_error = np.sqrt((1 + target**2) / nreal)
    assert np.mean(np.abs(cov - target) > 4 * standard_error) < 0.01


def test_seed():
    x, z = np.random.default_rng(0).random((2, 200))
    a = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1, seed=3)
    b = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1, seed=3)
    assert np.array_equal(a, b)
    np.random.seed(7)
    c = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1)
    np.random.seed(7)
    d = rf.get_random_field(x=x, z=z, theta_z=0.5, mu=0, sigma=1)
    assert np.array_equal(c, d)
