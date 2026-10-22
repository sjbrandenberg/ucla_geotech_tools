"""ipyconsol is a nonlinear implicit finite difference solver for one-dimensional consolidation of compressible soil with secondary compression.
ipyconsol is the Python version of the Javascript package iConsol.js published by Brandenberg (2017).
This code is an extension of Brandenberg (2017) because it the soil can be non-uniform and the loading can vary in time.
For spatially variable soil, material constants are entered as arrays, where the array length is equal to the number of nodes.
For time-varying loading, a loadfactor array is specified and has a length equal to the time array.

Brandenberg, S. J. (2017) "iConsol. js: JavaScript implicit finite-difference code for nonlinear consolidation and secondary compression." International Journal of Geomechanics, 17(6)

"""
import numpy as np
from scipy.linalg import solve_banded

LN10 = np.log(10.0)


def get_inputs(**kwargs):
    """Parse input arguments, initialize inputs, and return them as an array.

    There are two options for soil inputs (uniform and non-uniform) and two options for loading (constant and time-varying).

    soil property options:
        option 1: soil is uniform (i.e., material parameters are constants)
        option 2: soil is non-uniform (i.e., material parameters are arrays, with array lengths equal to number of nodes)

    load stage options:
        option A: constant load applied at time=0
        option B: time-varying load

    The options may be combined in any manner (e.g., [option 1, option A], [option 1, option B], [option 2, option A], [option 2, option B])
    Uniform soil has constant model coefficients, and either void ratio, OCR, or maximum past pressure is constant.

    Keyword Args:
        N: number of elements (required for option 1, scalar)
        H: thickness of soil layer (required for option 1, scalar)
        depth: array of depth values (required for option 2, array)
        Ntime: number of time steps (required for option A, scalar)
        tmax: maximum time value (required for option A, scalar)
        time: array of time values (required for option B, array)
        loadfactor: array of load factor values where stress increment = loadfactor*dsigv (required for option B, array).
        Cc: virgin compression index (required, scalar or array)
        Cr: recompression index (required, scalar or array)
        sigvref: vertical effective stress for a point on the normal consolidation line (required, scalar or array)
        esigvref: void ratio for point on normal consolidation line corresponding to esigvref (required, scalar or array)
        Gs: specific gravity of solids (required, scalar or array)
        kref: hydraulic conductivity for a point on the e-logk line (required, scalar or array)
        ekref: void ratio for point on the e-logk like corresponding to kref (required, scalar or array)
        Ck: coefficient of permeability variation defined as slope of e-logk relationship de/dlogk (required, scalar or array)
        Ca: coefficient of secondary compression (required, scalar or array)
        tref: time associated with normal consolidation line following Bjerrum's time-line concept (required, scalar or array)
        qo: initial vertical effective stress at top of soil (required, scalar)
        dsigv: change in total stress (required, scalar or array)
        ocrvoidratiotype: 0 OCR, 1 eo, 2 maximum past pressure (required, scalar or array)
        ocrvoidratio: value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype (required, scalar or array)
        ru: initial excess pore pressure ratio, prior to addition of dsigv (optional, scalar or array, default = 0)
        gammaw: unit weight of water (optional, default = 9.81)
        tol: convergence tolerance(optional, default = 1.0e-8)
        pa: atmospheric pressure (optional, default = 101.325)
        drainagetype: 0 = double-drained, 1 = single-drained through top, 2 = single-drained through bottom (optional, default = 0)
    """

    required = [
        "Cc",
        "Cr",
        "sigvref",
        "esigvref",
        "Gs",
        "kref",
        "ekref",
        "Ck",
        "Ca",
        "tref",
        "qo",
        "dsigv",
        "ocrvoidratiotype",
        "ocrvoidratio",
    ]
    for req in required:
        if req not in kwargs:
            print("You must specify " + req)
            return
    if np.ndim(kwargs["Cc"]) == 0:
        for req in required:
            if np.ndim(kwargs[req]) != 0:
                print("If Cc is a number, " + req + " must also be a number")
                return
        required1 = ["N", "H"]
        for req1 in required1:
            if req1 not in kwargs:
                print(
                    "you must specify "
                    + req1
                    + " when using input option 1, uniform soil"
                )
                return
        N = int(kwargs["N"])
        H = float(kwargs["H"])
        depth = np.linspace(0, H, num=N + 1, dtype=float, endpoint=True)
    elif isinstance(kwargs["Cc"], (list, tuple, np.ndarray)):
        for req in required:
            if req == "qo":
                continue
            if len(kwargs[req]) != len(kwargs["Cc"]):
                print(
                    "All vectors must have the same length using input option 2, layered soil. Cc has length = "
                    + str(len(kwargs["Cc"]))
                    + " while "
                    + req
                    + " has length = "
                    + str(len(kwargs[req]))
                )
                return
        if "depth" not in kwargs:
            print("You must specify depth if you specify Cc as a list or array.")
            return
        depth = np.asarray(kwargs["depth"], dtype=float)
        if len(depth) != len(kwargs["Cc"]):
            print("depth must have the same length as Cc")
            return
    else:
        print("Cc must be either a scalar, list, or np.ndarray")
        return

    # Scalars are broadcast to one value per node
    def nodal(value, dtype=float):
        return np.broadcast_to(np.asarray(value, dtype=dtype), depth.shape).copy()

    Cc, Cr, sigvref, esigvref, Gs, kref, ekref, Ck, Ca, tref = [
        nodal(kwargs[req]) for req in required[:10]
    ]
    qo = float(kwargs["qo"])
    dsigv = nodal(kwargs["dsigv"])
    ocrvoidratiotype = nodal(kwargs["ocrvoidratiotype"], dtype=np.int32)
    ocrvoidratio = nodal(kwargs["ocrvoidratio"])
    ru = nodal(kwargs.get("ru", 0.0))

    if "Ntime" in kwargs:
        if "tmax" not in kwargs:
            print("You must specify tmax if you specify Ntime")
            return
        else:
            time = np.logspace(
                np.log10(kwargs["tmax"]) - 5,
                np.log10(kwargs["tmax"]),
                num=kwargs["Ntime"],
                dtype=float,
            )
            loadfactor = np.full(kwargs["Ntime"], 1.0, dtype=float)
    elif "time" in kwargs:
        if "loadfactor" not in kwargs:
            print("You must specify loadfactor if you specify time")
            return
        if len(kwargs["loadfactor"]) != len(kwargs["time"]):
            print(
                "time and loadfactor must have the same length. Currently time has length = "
                + str(len(kwargs["time"]))
                + " and loadfactor has length = "
                + str(len(kwargs["loadfactor"]))
            )
            return
        time = np.asarray(kwargs["time"], dtype=float)
        loadfactor = np.asarray(kwargs["loadfactor"], dtype=float)
    else:
        print("You must specify either Ntime and tmax, or time and loadfactor")
        return

    gammaw = kwargs.get("gammaw", 9.81)
    tol = kwargs.get("tol", 1.0e-8)
    pa = kwargs.get("pa", 101.325)
    drainagetype = kwargs.get("drainagetype", 0)
    return [
        depth,
        Cc,
        Cr,
        sigvref,
        esigvref,
        Gs,
        kref,
        ekref,
        Ck,
        Ca,
        tref,
        qo,
        dsigv,
        ocrvoidratiotype,
        ocrvoidratio,
        ru,
        time,
        loadfactor,
        gammaw,
        tol,
        pa,
        drainagetype,
    ]


def _scalar_if_scalar(value, *inputs):
    """Return value as a float if every input is a scalar, so scalar calls return scalars"""
    if all(np.ndim(x) == 0 for x in inputs):
        return float(value)
    return value


def get_ktest(e, ekref, kref, Ck):
    """Return hydraulic conductivity value(s)

    Args:
        e = void ratio
        ekref = void ratio for any point on e-logk line
        kref = hydraulic conductivity for point on e-logk line corresponding to ekref
        Ck = coefficient of permeability variation (i.e., slope of e-logk line)
    """
    return kref * 10.0 ** ((np.asarray(e) - ekref) / Ck)


def _get_sigmap(elast, sigvlast, eref, sigvref, Cc, Cr):
    """Return maximum past pressure implied by void ratio and effective stress"""
    return 10.0 ** (
        (eref - elast + Cc * np.log10(sigvref) - Cr * np.log10(sigvlast)) / (Cc - Cr)
    )


def _creep(e, sigv, esigvref, sigvref, Cc, alpha):
    """Return exp((e - esigvref)/alpha + Cc/alpha*log10(sigv/sigvref)), which sets the rate of secondary
    compression following Bjerrum's time-line concept. Zero where alpha = 0 (no secondary compression)."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        value = np.exp(
            (e - esigvref) / alpha + Cc / alpha * np.log10(sigv / sigvref)
        )
    return np.where(alpha > 0, value, 0.0)


def get_avtest(elast, sigvlast, sigvtest, eref, sigvref, Cc, Cr):
    """Return compressibility av = -de/dsigv'

    Args (scalars or arrays):
        elast: void ratio at beginning of time step
        sigvlast: vertical effective stress at beginning of time step
        sigvtest: vertical effective stress at end of time step
        eref: void ratio for any point on normal consolidation line
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        Cc: virgin compression index
        Cr: recompression index
    """
    sigmap = _get_sigmap(elast, sigvlast, eref, sigvref, Cc, Cr)
    with np.errstate(divide="ignore", invalid="ignore"):
        virgin = (sigmap <= sigvlast) & (sigvtest > sigvlast)
        recompression = ~virgin & ((sigmap > sigvtest) | (sigvtest < sigvlast))
        De = np.where(
            virgin,
            Cc * np.log10(sigvtest / sigvlast),
            np.where(
                recompression,
                Cr * np.log10(sigvtest / sigvlast),
                Cr * np.log10(sigmap / sigvlast) + Cc * np.log10(sigvtest / sigmap),
            ),
        )
        unchanged = np.abs((sigvlast - sigvtest) / sigvtest) < 1.0e-5
        on_ncl = (np.abs((sigmap - sigvlast) / sigvlast) < 1.0e-5) & unchanged
        below_ncl = ~on_ncl & (sigmap > sigvtest) & unchanged
        av = np.where(
            on_ncl,
            Cc / sigvlast / LN10,
            np.where(below_ncl, Cr / sigvlast / LN10, De / (sigvtest - sigvlast)),
        )
    return _scalar_if_scalar(av, elast, sigvlast, sigvtest)


def get_etest(
    etest,
    elast,
    eref,
    alpha,
    tref,
    sigvtest,
    sigvlast,
    sigvref,
    Cc,
    avtest,
    utest,
    ulast,
    dt,
):
    """Returns trial value(s) of void ratio at the end of the current time step

    Newton-Raphson iterations are performed for all nodes with alpha > 0 at once, and each node stops
    iterating when its residual is below tolerance. A node that does not converge in 100 iterations
    returns -1.

    Args (scalars or arrays):
        etest: trial void ratio at the end of the time step
        elast: void ratio at the beginning of the time step
        eref: void ratio for any point on normal consolidation line
        alpha: secondary compression index (Ca/np.log(10.0))
        tref: time associated with normal consolidation line based on Bjerrum's time-line concept
        sigvtest: trial value of effective stress at the end of the time step
        sigvlast: value of effective stress at the beginning of the time step
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        Cc: virgin compression index
        avtest: trial compressibility at end of time step
        utest: trial excess pore pressure at end of time step
        ulast: excess pore ressure at beginning of time step
        dt: time step
    """
    etol = 1.0e-12
    inputs = np.broadcast_arrays(
        *[
            np.asarray(x, dtype=float)
            for x in (etest, elast, eref, alpha, tref, sigvtest, sigvlast, sigvref, Cc, avtest, utest, ulast)
        ]
    )
    shape = inputs[0].shape
    et, el, er, a, tr, st, sl, sr, cc, av, ut, ul = [x.ravel() for x in inputs]
    out = el + av * (ut - ul)
    creep = np.flatnonzero(a > 0)
    if creep.size:
        et, el, er, a, tr, st, sl, sr, cc, av, ut, ul = [
            x[creep] for x in (et, el, er, a, tr, st, sl, sr, cc, av, ut, ul)
        ]
        e = et.copy()
        last = np.exp((el - er) / a + cc / a * np.log10(sl / sr))

        def residual(e, k):
            return (
                av[k] * (ut[k] - ul[k])
                - 0.5
                * dt
                * a[k]
                / tr[k]
                * (np.exp((e[k] - er[k]) / a[k] + cc[k] / a[k] * np.log10(st[k] / sr[k])) + last[k])
                - e[k]
                + el[k]
            )

        res = residual(e, slice(None))
        active = np.flatnonzero(np.abs(res) > etol)
        for _ in range(100):
            if active.size == 0:
                break
            k = active
            dres_detest = (
                -0.5 * dt / tr[k] * np.exp((e[k] - er[k]) / a[k] + cc[k] / a[k] * np.log10(st[k] / sr[k]))
                - 1.0
            )
            e[k] = e[k] - res[k] / dres_detest
            res[k] = residual(e, k)
            active = k[np.abs(res[k]) > etol]
        e[active] = -1.0
        out[creep] = e
    return _scalar_if_scalar(out.reshape(shape), etest, elast, avtest, utest)


def _flux(u, k, dz2, i, j):
    """Return the finite difference flux terms of the residual, per dt/gammaw, at nodes i with neighbors
    j on one side and an impermeable boundary on the other (used for the undrained boundary node)"""
    return 0.5 * k[i] * (u[j] - u[i]) / dz2 + 0.125 * (u[i] - u[j]) * (k[i] - k[j]) / dz2


def get_residual(
    ulast,
    utest,
    elast,
    etest,
    klast,
    ktest,
    avtest,
    zlast,
    ztest,
    sigvlast,
    sigvtest,
    gammaw,
    Ca,
    Cc,
    sigvref,
    esigvref,
    dt,
    tref,
    drainagetype,
    N,
):
    """Returns an array of residuals

    Args:
        ulast: array of excess pore pressures at beginning of time step
        utest: array of trial excess pore pressures at end of time step
        elast: array of void ratios at the beginning of the time step
        etest: array of trial void ratios at the end of the time step
        klast: array of hydraulic conductivities at beginning of time step
        ktest: array of trial hydraulic conductivities at end of time step
        avtest: array of trial coefficients of compressibility at end of time step
        zlast: array of node depths at beginning of time step
        ztest: array of trial node depths at end of time step
        sigvlast: array of values of effective stress at the beginning of the time step
        sigvtest: array of trial values of effective stress at the end of the time step
        gammaw: unit weight of water
        Ca: array of secondary compression coefficients
        Cc: array of virgin compression indices
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        esigvref: void ratio for any point on normal consolidation line
        dt: time step
        tref: array of times associated with normal consolidation line based on Bjerrum's time-line concept
        drainagetype: 0 = double-drained, 1 = single drained through top, 2 = single-drained through bottom
        N: number of elements
    """
    alpha = Ca / LN10
    s, m, p = slice(1, N), slice(0, N - 1), slice(2, N + 1)
    dzlast2 = (0.5 * (zlast[p] - zlast[m])) ** 2
    dztest2 = (0.5 * (ztest[p] - ztest[m])) ** 2
    residual = -avtest / (1.0 + elast) * (utest - ulast)
    residual += (
        0.5
        * alpha
        * dt
        / tref
        / (1.0 + elast)
        * (
            _creep(etest, sigvtest, esigvref, sigvref, Cc, alpha)
            + _creep(elast, sigvlast, esigvref, sigvref, Cc, alpha)
        )
    )
    residual[s] += (
        0.5 * dt / gammaw * ktest[s] * (utest[p] - 2.0 * utest[s] + utest[m]) / dztest2
        + 0.5 * dt / gammaw * klast[s] * (ulast[p] - 2.0 * ulast[s] + ulast[m]) / dzlast2
        + 0.125 * dt / gammaw * (utest[p] - utest[m]) * (ktest[p] - ktest[m]) / dztest2
        + 0.125 * dt / gammaw * (ulast[p] - ulast[m]) * (klast[p] - klast[m]) / dzlast2
    )
    if drainagetype == 1:
        residual[N] += (
            dt
            / gammaw
            * (
                _flux(utest, ktest, (ztest[N] - ztest[N - 1]) ** 2, N, N - 1)
                + _flux(ulast, klast, (zlast[N] - zlast[N - 1]) ** 2, N, N - 1)
            )
        )
    if drainagetype == 2:
        residual[0] += (
            dt
            / gammaw
            * (
                _flux(utest, ktest, (ztest[1] - ztest[0]) ** 2, 0, 1)
                + _flux(ulast, klast, (zlast[1] - zlast[0]) ** 2, 0, 1)
            )
        )
    # Drained boundaries have prescribed pore pressure
    if drainagetype in (0, 1):
        residual[0] = 0.0
    if drainagetype in (0, 2):
        residual[N] = 0.0
    return residual


def _solve_tridiagonal(lower, diag, upper, rhs, drainagetype, N):
    """Solve the tridiagonal system for the nodes whose pore pressure is unknown.

    lower, diag, and upper are the coefficients of u[i-1], u[i], and u[i+1] in the equation for node i.
    Drained boundary nodes are excluded, and zero is returned at those nodes."""
    lo, hi = {0: (1, N), 1: (1, N + 1), 2: (0, N)}[drainagetype]
    ab = np.zeros((3, hi - lo))
    ab[0, 1:] = upper[lo : hi - 1]
    ab[1] = diag[lo:hi]
    ab[2, :-1] = lower[lo + 1 : hi]
    x = np.zeros(N + 1)
    x[lo:hi] = solve_banded((1, 1), ab, rhs[lo:hi])
    return x


def get_utest2(
    N,
    klast,
    ktest,
    zlast,
    ztest,
    avtest,
    elast,
    etest,
    ulast,
    utest,
    sigvlast,
    sigvtest,
    Res,
    dt,
    Cc,
    Cr,
    Ca,
    gammaw,
    ekref,
    kref,
    Ck,
    tref,
    esigvref,
    sigvref,
    dsigv,
    drainagetype,
    pa,
):
    """Returns an array of pore pressures from one Newton-Raphson iteration on the residuals

    ktest, ztest, avtest, sigvtest, and pa are not used, and are retained for backward compatibility.

    Args:
        N: number of elements
        klast: array of hydraulic conductivities at beginning of time step
        ktest: array of trial hydraulic conductivities at end of time step
        zlast: array of node depths at beginning of time step
        ztest: array of trial node depths at end of time step
        avtest: array of trial coefficients of compressibility at end of time step
        elast: array of void ratios at the beginning of the time step
        etest: array of trial void ratios at the end of the time step
        ulast: array of excess pore pressures at beginning of time step
        utest: array of trial excess pore pressures at end of time step
        sigvlast: array of values of effective stress at the beginning of the time step
        sigvtest: array of trial values of effective stress at the end of the time step
        Res: array of residuals
        dt: time step
        Cc: array of virgin compression indices
        Cr: array of recompression indices
        Ca: array of secondary compression coefficients
        gammaw: unit weight of water
        kref = hydraulic conductivity for any point on e-logk line
        ekref = void ratio for point on e-logk line corresponding to kref
        Ck = coefficient of permeability variation (i.e., slope of e-logk line)
        tref: array of times associated with normal consolidation line based on Bjerrum's time-line concept
        esigvref: void ratio for any point on normal consolidation line
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        dsigv: change in vertical total stress
        drainagetype: 0 = double-drained, 1 = single drained through top, 2 = single-drained through bottom
        pa: atmospheric pressure
    """
    # The residual at node i depends only on pore pressures at nodes i-1, i, and i+1, so the Jacobian is
    # tridiagonal. It is found by finite differences, perturbing every third node at once.
    def residual(u):
        sigv, av, e, k, z = _get_trial_state(
            u, ulast, sigvlast, elast, etest, Cc, Cr, Ca, esigvref, sigvref, tref, ekref, kref, Ck, zlast, zlast[-1], dt
        )
        return get_residual(
            ulast, u, elast, e, klast, k, av, zlast, z, sigvlast, sigv, gammaw, Ca, Cc, sigvref, esigvref, dt, tref,
            drainagetype, N,
        )

    h = 1.0e-7 * (np.abs(utest) + np.max(np.abs(dsigv)) + 1.0)
    lower = np.zeros(N + 1)
    diag = np.zeros(N + 1)
    upper = np.zeros(N + 1)
    with np.errstate(all="ignore"):
        for color in range(3):
            j = np.arange(color, N + 1, 3)
            uplus = utest.copy()
            uplus[j] += h[j]
            dRes = residual(uplus) - Res
            diag[j] = dRes[j] / h[j]
            below = j[j < N]
            lower[below + 1] = dRes[below + 1] / h[below]
            above = j[j > 0]
            upper[above - 1] = dRes[above - 1] / h[above]
    return utest - _solve_tridiagonal(lower, diag, upper, Res, drainagetype, N)


def get_utest(
    klast,
    ktest,
    avtest,
    zlast,
    ztest,
    ulast,
    sigvlast,
    sigvtest,
    elast,
    etest,
    Ca,
    Cc,
    sigvref,
    esigvref,
    dt,
    tref,
    drainagetype,
    N,
    gammaw,
):
    """Returns an initial guess of array of pore pressures based on current state of the model

    Args:
        klast: array of hydraulic conductivities at beginning of time step
        ktest: array of trial hydraulic conductivities at end of time step
        avtest: array of trial coefficients of compressibility at end of time step
        zlast: array of node depths at beginning of time step
        ztest: array of trial node depths at end of time step
        ulast: array of excess pore pressures at beginning of time step
        sigvlast: array of values of effective stress at the beginning of the time step
        sigvtest: array of trial values of effective stress at the end of the time step
        elast: array of void ratios at the beginning of the time step
        etest: array of trial void ratios at the end of the time step
        Ca: array of secondary compression coefficients
        Cc: array of virgin compression indices
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        esigvref: void ratio for any point on normal consolidation line
        dt: time step
        tref: array of times associated with normal consolidation line based on Bjerrum's time-line concept
        drainagetype: 0 = double-drained, 1 = single drained through top, 2 = single-drained through bottom
        N: number of elements
        gammaw: unit weight of water
    """
    alpha = Ca / LN10
    k, kl, ul = ktest, klast, ulast
    lower = np.zeros(N + 1)
    upper = np.zeros(N + 1)
    diag = avtest / (1.0 + elast)
    rhs = avtest * ulast / (1.0 + elast) + 0.5 * alpha * dt / tref / (1.0 + elast) * (
        _creep(elast, sigvlast, esigvref, sigvref, Cc, alpha)
        + _creep(etest, sigvtest, esigvref, sigvref, Cc, alpha)
    )

    # Interior nodes
    s, m, p = slice(1, N), slice(0, N - 1), slice(2, N + 1)
    g = 0.5 * dt / gammaw / (0.5 * (ztest[p] - ztest[m])) ** 2
    glast = 0.5 * dt / gammaw / (0.5 * (zlast[p] - zlast[m])) ** 2
    lower[s] = -g * (k[s] + 0.25 * (k[m] - k[p]))
    diag[s] += 2.0 * g * k[s]
    upper[s] = -g * (k[s] - 0.25 * (k[m] - k[p]))
    rhs[s] += glast * kl[s] * (ul[m] - 2.0 * ul[s] + ul[p]) + 0.25 * glast * (ul[p] - ul[m]) * (kl[p] - kl[m])

    # Undrained boundary node
    if drainagetype in (1, 2):
        i, j = (N, N - 1) if drainagetype == 1 else (0, 1)
        g = 0.5 * dt / gammaw / (ztest[i] - ztest[j]) ** 2
        glast = 0.5 * dt / gammaw / (zlast[i] - zlast[j]) ** 2
        diag[i] += g * k[i]
        neighbor = -g * (k[i] + 0.25 * (k[j] - k[i]))
        if drainagetype == 1:
            lower[i] = neighbor
        else:
            upper[i] = neighbor
        rhs[i] += glast * kl[i] * (ul[j] - ul[i]) + 0.25 * glast * (ul[i] - ul[j]) * (kl[i] - kl[j])

    return _solve_tridiagonal(lower, diag, upper, rhs, drainagetype, N)


def Max(v):
    """Returns the maximum of the absolute value of an array

    Args:
        v: array of values
    """
    return np.max(np.abs(v))


def get_initial(
    depth,
    Cc,
    Cr,
    sigvref,
    esigvref,
    Gs,
    kref,
    ekref,
    Ck,
    Ca,
    tref,
    qo,
    dsigv,
    ocrvoidratiotype,
    ocrvoidratio,
    ru,
    time,
    loadfactor,
    gammaw,
    tol,
    pa,
    drainagetype,
):
    """Returns Python dictionary containing initial values of void ratio and vertial effective stress, and final vertical effective stress

    The initial state is integrated from the top of the layer downward, so this calculation proceeds node by node.

    Args:
        depth: array of nodal depths
        Cc: array of virgin compression indices
        Cr: array of recompression indices
        sigvref: vertical effective stress for point on normal consolidation line corresponding to eref
        esigvref: void ratio for any point on normal consolidation line
        Gs: array of specific gravity values
        kref = array of hydraulic conductivity for any point on e-logk line
        ekref = array of void ratios for point on e-logk line corresponding to kref
        Ck = array of coefficients of permeability variation (i.e., slope of e-logk line)
        Ca = array of secondary compression indices
        tref: array of times associated with normal consolidation line based on Bjerrum's time-line concept
        qo: initial vertical effective stress at surface of compressible layer
        dsigv: change in vertical total stress
        ocrvoidratiotype = 0 constant OCR, 1 constant eo, 2 constant maximum past pressure
        ocrvoidratio = value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype
        ru = array of initial excess pore pressure ratios, prior to addition of dsigv
        time = array of time values
        loadfactor = corresponding load factor values. Stress increment is load factor multiplied by dsigv.
        gammaw: unit weight of water
        tol: convergence tolerance
        pa: atmospheric pressure
        drainagetype: 0 = double-drained, 1 = single drained through top, 2 = single-drained through bottom
    """
    sigvo = np.zeros(depth.size, dtype=np.float64)
    eo = np.zeros(depth.size, dtype=np.float64)

    def void_ratio(i, sigv):
        if ocrvoidratiotype[i] == 0:
            return (
                esigvref[i]
                - Cc[i] * np.log10(ocrvoidratio[i] * sigv / sigvref[i])
                + Cr[i] * np.log10(ocrvoidratio[i])
            )
        if ocrvoidratiotype[i] == 1:
            return ocrvoidratio[i]
        return (
            esigvref[i]
            - Cc[i] * np.log10(ocrvoidratio[i] / sigvref[i])
            + Cr[i] * np.log10(ocrvoidratio[i] / sigv)
        )

    for i in range(depth.size):
        if ocrvoidratiotype[i] not in (0, 1, 2):
            print("Invalid ocrvoidratiotype: " + str(ocrvoidratiotype[i]))
            return
    sigvo[0] = qo
    eo[0] = void_ratio(0, qo)
    for i in range(1, depth.size):
        dz = depth[i] - depth[i - 1]
        if ocrvoidratiotype[i] == 1:
            eo[i] = ocrvoidratio[i]
            sigvo[i] = sigvo[i - 1] + (Gs[i] - 1.0) * dz * gammaw / (
                1.0 + 0.5 * (eo[i] + eo[i - 1])
            )
            continue
        # Fixed-point iteration on effective stress, using the average void ratio over the element
        sigvotest = sigvo[i - 1] + (Gs[i] - 1.0) / (1.0 + eo[i - 1]) * gammaw * dz
        eotest = void_ratio(i, sigvotest)
        res = (
            sigvotest
            - sigvo[i - 1]
            - (Gs[i] - 1.0) * gammaw * dz / (1.0 + 0.5 * (eotest + eo[i - 1]))
        )
        while np.fabs(res) > 1.0e-10 * pa:
            sigvotest = sigvo[i - 1] + (Gs[i] - 1.0) * gammaw * dz / (
                1.0 + 0.5 * (eotest + eo[i - 1])
            )
            eotest = void_ratio(i, sigvotest)
            res = (
                sigvotest
                - sigvo[i - 1]
                - (Gs[i] - 1.0) * gammaw * dz / (1.0 + 0.5 * (eotest + eo[i - 1]))
            )
        sigvo[i] = sigvotest
        eo[i] = eotest
    sigvo2 = sigvo * (1.0 - ru)
    sigvf = sigvo + dsigv * loadfactor[0]
    return {"eo": eo, "sigvo": sigvo, "sigvo2": sigvo2, "sigvf": sigvf}


def _get_trial_state(
    utest, ulast, sigvlast, elast, etest, Cc, Cr, Ca, esigvref, sigvref, tref, ekref, kref, Ck, zlast, H, dt
):
    """Return effective stress, compressibility, void ratio, hydraulic conductivity, and node depths
    at the end of the time step for trial pore pressures utest"""
    sigvtest = sigvlast + ulast - utest
    avtest = get_avtest(elast, sigvlast, sigvtest, esigvref, sigvref, Cc, Cr)
    etest = get_etest(
        etest, elast, esigvref, Ca / LN10, tref, sigvtest, sigvlast, sigvref, Cc, avtest, utest, ulast, dt
    )
    epsvtest = (elast - etest) / (1.0 + elast)
    ktest = get_ktest(etest, ekref, kref, Ck)
    # Node depths follow from compression of each element, with the bottom of the layer fixed
    ztest = np.empty_like(zlast)
    ztest[-1] = H
    ztest[:-1] = H - np.cumsum((np.diff(zlast) * (1.0 - epsvtest[:-1]))[::-1])[::-1]
    return sigvtest, avtest, etest, ktest, ztest


def compute(**kwargs):
    """Main function for solving nonlinear consolidation problem. Returns a Python dictionary containing depths, excess pore pressures,
    vertical effective stresses, and void ratios at each time and depth

    There are two options for soil inputs (uniform and non-uniform) and two options for loading (constant and time-varying).

    soil property options:
        option 1: soil is uniform (i.e., material parameters are constants)
        option 2: soil is non-uniform (i.e., material parameters are arrays, with array lengths equal to number of nodes)

    load stage options:
        option A: constant load applied at time=0
        option B: time-varying load

    The options may be combined in any manner (e.g., [option 1, option A], [option 1, option B], [option 2, option A], [option 2, option B])
    Uniform soil has constant model coefficients, and either void ratio, OCR, or maximum past pressure is constant.

    Keyword Args:
        N: number of elements (required for option 1, scalar)
        H: thickness of soil layer (required for option 1, scalar)
        depth: array of depth values (required for option 2, array)
        Ntime: number of time steps (required for option A, scalar)
        tmax: maximum time value (required for option A, scalar)
        time: array of time values (required for option B, array)
        loadfactor: array of load factor values where stress increment = loadfactor*dsigv (required for option B, array).
        Cc: virgin compression index (required, scalar or array)
        Cr: recompression index (required, scalar or array)
        sigvref: vertical effective stress for a point on the normal consolidation line (required, scalar or array)
        esigvref: void ratio for point on normal consolidation line corresponding to esigvref (required, scalar or array)
        Gs: specific gravity of solids (required, scalar or array)
        kref: hydraulic conductivity for a point on the e-logk line (required, scalar or array)
        ekref: void ratio for point on the e-logk like corresponding to kref (required, scalar or array)
        Ck: coefficient of permeability variation defined as slope of e-logk relationship de/dlogk (required, scalar or array)
        Ca: coefficient of secondary compression (required, scalar or array)
        tref: time associated with normal consolidation line following Bjerrum's time-line concept (required, scalar or array)
        qo: initial vertical effective stress at top of soil (required, scalar)
        dsigv: change in total stress (required, scalar or array)
        ocrvoidratiotype: 0 OCR, 1 eo, 2 maximum past pressure (required, scalar or array)
        ocrvoidratio: value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype (required, scalar or array)
        ru: initial excess pore pressure ratio, prior to addition of dsigv (optional, scalar or array, default = 0)
        gammaw: unit weight of water (optional, default = 9.81)
        tol: convergence tolerance(optional, default = 1.0e-8)
        pa: atmospheric pressure (optional, default = 101.325)
        drainagetype: 0 = double-drained, 1 = single-drained through top, 2 = single-drained through bottom (optional, default = 0)
    """
    inputs = get_inputs(**kwargs)
    if inputs is None:
        return
    (
        depth,
        Cc,
        Cr,
        sigvref,
        esigvref,
        Gs,
        kref,
        ekref,
        Ck,
        Ca,
        tref,
        qo,
        dsigv,
        ocrvoidratiotype,
        ocrvoidratio,
        ru,
        time,
        loadfactor,
        gammaw,
        tol,
        pa,
        drainagetype,
    ) = inputs

    N = depth.size - 1
    Ntime = time.size
    H = depth[N]

    u = np.zeros((N + 1, Ntime), dtype=float)
    sigv = np.zeros((N + 1, Ntime), dtype=float)
    z = np.zeros((N + 1, Ntime), dtype=float)
    e = np.zeros((N + 1, Ntime), dtype=float)

    initial = get_initial(
        depth,
        Cc,
        Cr,
        sigvref,
        esigvref,
        Gs,
        kref,
        ekref,
        Ck,
        Ca,
        tref,
        qo,
        dsigv,
        ocrvoidratiotype,
        ocrvoidratio,
        ru,
        time,
        loadfactor,
        gammaw,
        tol,
        pa,
        drainagetype,
    )
    if initial is None:
        return

    dsigvsublast = dsigv * loadfactor[0]
    ulast = initial["sigvf"] - initial["sigvo2"]
    sigvlast = initial["sigvf"] - ulast
    elast = initial["eo"].copy()
    klast = get_ktest(elast, ekref, kref, Ck)
    zlast = depth.copy()
    material = dict(
        Cc=Cc, Cr=Cr, Ca=Ca, esigvref=esigvref, sigvref=sigvref, tref=tref, ekref=ekref, kref=kref, Ck=Ck, H=H
    )

    for j in range(Ntime):
        dsigvsub = dsigv * loadfactor[j]
        ulast = ulast + dsigvsub - dsigvsublast
        dsigvsublast = dsigvsub
        dt = time[j] if j == 0 else time[j] - time[j - 1]
        ktest = klast.copy()
        etest = elast.copy()
        avtest = get_avtest(elast, sigvlast, sigvlast, esigvref, sigvref, Cc, Cr)
        utest = get_utest(
            klast,
            ktest,
            avtest,
            zlast,
            zlast,
            ulast,
            sigvlast,
            sigvlast,
            elast,
            etest,
            Ca,
            Cc,
            sigvref,
            esigvref,
            dt,
            tref,
            drainagetype,
            N,
            gammaw,
        )

        def trial(utest, etest):
            """Return the state at the end of the time step for trial pore pressures, and its residuals"""
            sigvtest, avtest, etest, ktest, ztest = _get_trial_state(
                utest, ulast, sigvlast, elast, etest, zlast=zlast, dt=dt, **material
            )
            Res = get_residual(
                ulast,
                utest,
                elast,
                etest,
                klast,
                ktest,
                avtest,
                zlast,
                ztest,
                sigvlast,
                sigvtest,
                gammaw,
                Ca,
                Cc,
                sigvref,
                esigvref,
                dt,
                tref,
                drainagetype,
                N,
            )
            return (sigvtest, avtest, etest, ktest, ztest), Res

        with np.errstate(all="ignore"):
            state, Res = trial(utest, etest)
            if not np.all(np.isfinite(Res)):
                # The predictor overshot (e.g., negative effective stress), so start from the previous pore pressures
                utest = ulast.copy()
                state, Res = trial(utest, etest)
        I = 0
        # NaN residuals compare False with tol, so check them explicitly
        while not (np.all(np.isfinite(Res)) and Max(Res) <= tol):
            if I == 100 or not np.all(np.isfinite(Res)):
                print(
                    "ERROR! Increment "
                    + str(j)
                    + " did not converge in "
                    + str(I)
                    + " iterations."
                )
                return {"z": z, "u": u, "sigv": sigv, "e": e}
            sigvtest, avtest, etest, ktest, ztest = state
            step = get_utest2(
                N,
                klast,
                ktest,
                zlast,
                ztest,
                avtest,
                elast,
                etest,
                ulast,
                utest,
                sigvlast,
                sigvtest,
                Res,
                dt,
                Cc,
                Cr,
                Ca,
                gammaw,
                ekref,
                kref,
                Ck,
                tref,
                esigvref,
                sigvref,
                dsigvsub,
                drainagetype,
                pa,
            ) - utest
            # Backtracking line search: halve the Newton step until the state is valid and the residual decreases
            lam = 1.0
            with np.errstate(all="ignore"):
                for _ in range(40):
                    candidate = utest + lam * step
                    candidate_state, candidate_Res = trial(candidate, etest)
                    if np.all(np.isfinite(candidate_Res)) and Max(candidate_Res) < Max(Res):
                        break
                    lam *= 0.5
            utest, state, Res = candidate, candidate_state, candidate_Res
            I += 1
        sigvtest, avtest, etest, ktest, ztest = state
        u[:, j] = ulast = utest
        sigv[:, j] = sigvlast = sigvtest
        e[:, j] = elast = etest
        z[:, j] = zlast = ztest
        klast = ktest
    return {"z": z, "u": u, "sigv": sigv, "e": e}
