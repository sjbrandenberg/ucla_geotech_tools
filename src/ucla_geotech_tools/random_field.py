"""random_field is a tool for efficiently creating realizations of two-dimensional random fields by combining Cholesky decomposition of a covariance matrix with kriging at intermediate nodes.

Details of the method are described by Yang and Brandenberg (2022).

Yang, Y., Wang, P., and Brandenberg, S. J. (2022) "An algorithm for generating spatially correlated random fields using Cholesky decomposition and ordinary kriging." Computers and Geotechnics, 147

"""
import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.spatial import cKDTree

# Added to the diagonal of the correlation matrix at Cholesky nodes so that it is positive definite and
# well conditioned. The same matrix is used for the Cholesky decomposition and for kriging.
NUGGET = 1.0e-4

# Default maximum distance from a kriging node to the nearest Cholesky node, as a multiple of the correlation
# length. Kriging reproduces a smooth (Gaussian) field closely from widely spaced Cholesky nodes, but a rough
# (exponential) field requires closely spaced Cholesky nodes.
DHCHOL_FACTOR = {0: 0.15, 1: 0.07}

# Maximum number of entries in a block of the correlation matrix between Cholesky and kriging nodes
_BLOCK_SIZE = 2**22


def get_C(x, z, theta, model):
    """Returns a correlation matrix for an input array containing x and z coordinates, a correlation length, and the type of correlation model (0 = Gaussian, 1 = exponential)

    Args:
        x: array containing x coordinates
        z: array containing z coordinates
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    dist = np.sqrt((x - x[:, np.newaxis]) ** 2.0 + (z - z[:, np.newaxis]) ** 2.0)
    if model == 0:
        C = np.exp(-np.pi * (dist / theta) ** 2.0)
    else:
        C = np.exp(-2.0 * dist / theta)

    return C


def get_C_chol_krig(xchol, zchol, xkrig, zkrig, theta, model):
    """Returns correlation matrix between the cholesky decomposition nodes and the kriging nodes, with one row per kriging node.

    Args:
        xchol: array containing x coordinates of cholesky nodes
        zchol: array containing z coordinates of cholesky nodes
        xkrig: array containing x coordinates of kriging nodes
        zkrig: array containing z coordinates of kriging nodes
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    dist = np.sqrt((xchol - xkrig[:, np.newaxis]) ** 2.0 + (zchol - zkrig[:, np.newaxis]) ** 2.0)
    if model == 0:
        C = np.exp(-np.pi * (dist / theta) ** 2.0)
    else:
        C = np.exp(-2.0 * dist / theta)

    return C


def get_chol_krig_nodes(x, z, dist_max):
    """Divides nodes into a set for cholesky decomposition and a set for kriging.

    Nodes are visited in order. A node becomes a kriging node if it is within dist_max of a cholesky node,
    and otherwise becomes a cholesky node. Every node is assigned to exactly one set, and every kriging node
    is within dist_max of a cholesky node.

    Args:
    x: x coordinates of all nodes
    z: z coordinates of all nodes
    dist_max: maximum distance from a kriging node to the nearest cholesky node
    """
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    tree = cKDTree(np.column_stack((x, z)))
    assigned = np.zeros(len(x), dtype=bool)
    chol = []
    for i in range(len(x)):
        if assigned[i]:
            continue
        chol.append(i)
        assigned[tree.query_ball_point((x[i], z[i]), dist_max, p=2.0, eps=0.0)] = True
    ind_chol = np.asarray(chol, dtype=int)
    is_chol = np.zeros(len(x), dtype=bool)
    is_chol[ind_chol] = True
    ind_krig = np.flatnonzero(~is_chol)
    return (ind_chol, x[ind_chol], z[ind_chol], ind_krig, x[ind_krig], z[ind_krig])


def _standard_normal(size, rng):
    """Draw standard normal variates from rng, or from numpy's global generator if rng is None"""
    if rng is None:
        return np.random.normal(loc=0, scale=1, size=size)
    return rng.standard_normal(size)


def get_chol_field(C, mu, sigma, rng=None):
    """Function that uses cholesky decomposition to return a random realization at the cholesky nodes

    Args:
        C: correlation matrix at cholesky nodes
        mu: mean of random field
        sigma: standard deviation of random field
        rng: numpy random Generator (optional; numpy's global generator is used if None)
    """
    Nchol = len(C)
    rv = _standard_normal(Nchol, rng)
    L_chol = np.linalg.cholesky(C + NUGGET * np.identity(Nchol))
    chol_field = mu + sigma * L_chol @ rv

    return chol_field


def _kriging_weights(L, chol_field):
    """Return the generalized least squares mean and the vector alpha such that the ordinary kriging estimate at
    a node with correlation vector c to the cholesky nodes is mean + c @ alpha.

    Ordinary kriging is evaluated as simple kriging about the generalized least squares estimate of the mean,
    which gives the same estimates as solving the ordinary kriging system, using the cholesky factor L of the
    correlation matrix at cholesky nodes.
    """
    ones = solve_triangular(L, np.ones(len(L)), lower=True)
    mean = ones @ solve_triangular(L, chol_field, lower=True) / (ones @ ones)
    alpha = cho_solve((L, True), chol_field - mean)
    return mean, alpha


def _correlation(dist, theta, model):
    if model == 0:
        return np.exp(-np.pi * (dist / theta) ** 2.0)
    return np.exp(-2.0 * dist / theta)


def _kriging_variance(xchol, zchol, xkrig, zkrig, theta, model, tree, neighbors=32):
    """Return ordinary kriging variances (for unit variance) at kriging nodes, computed from the nearest
    cholesky nodes. Kriging weights decay rapidly beyond the nearest cholesky nodes, so this closely
    approximates the variance from all cholesky nodes at a small fraction of the cost."""
    m = min(neighbors, len(xchol))
    _, nb = tree.query(np.column_stack((xkrig, zkrig)), k=m)
    nb = nb.reshape(len(xkrig), m)
    xn, zn = xchol[nb], zchol[nb]
    K = np.ones((len(xkrig), m + 1, m + 1))
    K[:, :m, :m] = _correlation(
        np.hypot(xn[:, :, None] - xn[:, None, :], zn[:, :, None] - zn[:, None, :]), theta, model
    ) + NUGGET * np.identity(m)
    K[:, m, m] = 0.0
    c = np.ones((len(xkrig), m + 1))
    c[:, :m] = _correlation(np.hypot(xn - xkrig[:, None], zn - zkrig[:, None]), theta, model)
    w = np.linalg.solve(K, c[:, :, None])[:, :, 0]
    return np.clip(1.0 - np.sum(w * c, axis=1), 0.0, None)


def get_krig_field(C, C_chol_krig, chol_field):
    """Returns values of random field at kriging nodes, which are ordinary kriging estimates from the cholesky nodes

    Args:
        C: correlation matrix at cholesky nodes
        C_chol_krig: correlation matrix between cholesky nodes and kriging nodes (one row per kriging node)
        chol_field: realization of random field at cholesky nodes
    """
    L = np.linalg.cholesky(C + NUGGET * np.identity(len(C)))
    mean, alpha = _kriging_weights(L, chol_field)
    return mean + C_chol_krig @ alpha


def get_random_field(**kwargs):
    """Function that returns composite random field from cholesky decomposition combined with kriging.

    Kriging estimates are smoother than the random field, so an independent normal residual with the
    kriging variance is added at each kriging node. The variance of the field is then sigma^2 at every node.
    The correlation between pairs of kriging nodes is slightly underestimated when the cholesky nodes are
    widely spaced relative to the correlation length, and the default Dhchol limits that error to about 0.03.

    Keyword Args:
        x: x coordinates of all nodes (required)
        z: z coordinates of all nodes (required)
        theta_z: correlation length in z-direction (required)
        mu: mean of random field (required)
        sigma: standard deviation of random field (required)
        model: correlation model, 0 = Gaussian 1 = exponential (default 0)
        theta_x: correlation length in x-direction (default theta_z)
        Dhchol: maximum distance from a kriging node to the nearest cholesky node, in z-direction units after
            scaling x by theta_z/theta_x (default 0.15*theta_z for Gaussian, 0.07*theta_z for exponential)
        seed: seed or numpy random Generator for reproducible fields (optional; numpy's global generator is
            used if not specified)
    """
    required = ["x", "z", "theta_z", "mu", "sigma"]
    for r in required:
        if r not in kwargs.keys():
            print(r + ' must be included like this "' + r + '= value"')
            return
    x = np.asarray(kwargs["x"], dtype=float)
    z = np.asarray(kwargs["z"], dtype=float)
    theta_z = kwargs["theta_z"]
    mu = kwargs["mu"]
    sigma = kwargs["sigma"]
    model = kwargs.get("model", 0)
    theta_x = kwargs.get("theta_x", theta_z)
    Dhchol = kwargs.get("Dhchol", DHCHOL_FACTOR[0 if model == 0 else 1] * theta_z)
    rng = np.random.default_rng(kwargs["seed"]) if "seed" in kwargs else None

    x = x * theta_z / theta_x  # perform coordinate transformation
    chol_ind, xvchol, zvchol, krig_ind, xvkrig, zvkrig = get_chol_krig_nodes(x, z, Dhchol)
    C = get_C(xvchol, zvchol, theta_z, model)
    L = np.linalg.cholesky(C + NUGGET * np.identity(len(C)))
    rf = np.zeros(len(x))
    rf[chol_ind] = L @ _standard_normal(len(chol_ind), rng)
    mean, alpha = _kriging_weights(L, rf[chol_ind])
    tree = cKDTree(np.column_stack((xvchol, zvchol)))
    # Krige in blocks of nodes to bound memory
    block = max(1, min(_BLOCK_SIZE // max(1, len(chol_ind)), 4096))
    for start in range(0, len(krig_ind), block):
        stop = start + block
        C_chol_krig = get_C_chol_krig(xvchol, zvchol, xvkrig[start:stop], zvkrig[start:stop], theta_z, model)
        variance = _kriging_variance(xvchol, zvchol, xvkrig[start:stop], zvkrig[start:stop], theta_z, model, tree)
        residual = np.sqrt(variance) * _standard_normal(len(variance), rng)
        rf[krig_ind[start:stop]] = mean + C_chol_krig @ alpha + residual

    return mu + sigma * rf
