import numpy as np
from scipy.signal import lfilter, lfiltic
from scipy.spatial import ConvexHull, QhullError


def get_ngawest2_T():
    """Return numpy array containing spectral periods utilized in NGAWest2 project"""
    return np.asarray(
        [
            0.01,
            0.02,
            0.022,
            0.025,
            0.029,
            0.030,
            0.032,
            0.035,
            0.036,
            0.040,
            0.042,
            0.044,
            0.045,
            0.046,
            0.048,
            0.050,
            0.055,
            0.060,
            0.065,
            0.067,
            0.070,
            0.075,
            0.080,
            0.085,
            0.090,
            0.095,
            0.100,
            0.11,
            0.12,
            0.13,
            0.133,
            0.14,
            0.15,
            0.16,
            0.17,
            0.18,
            0.19,
            0.20,
            0.22,
            0.24,
            0.25,
            0.26,
            0.28,
            0.29,
            0.30,
            0.32,
            0.34,
            0.35,
            0.36,
            0.38,
            0.40,
            0.42,
            0.44,
            0.45,
            0.46,
            0.48,
            0.50,
            0.55,
            0.60,
            0.65,
            0.667,
            0.70,
            0.75,
            0.80,
            0.85,
            0.90,
            0.95,
            1.0,
            1.1,
            1.2,
            1.3,
            1.4,
            1.5,
            1.6,
            1.7,
            1.8,
            1.9,
            2.0,
            2.2,
            2.4,
            2.5,
            2.6,
            2.8,
            3.0,
            3.2,
            3.4,
            3.5,
            3.6,
            3.8,
            4.0,
            4.2,
            4.4,
            4.6,
            4.8,
            5.0,
            5.5,
            6.0,
            6.5,
            7.0,
            7.5,
            8.0,
            8.5,
            9.0,
            9.5,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            20.0,
        ]
    )


def next_fast_len(target):
    """Return a number of datapoints for efficient calculation of fast Fourier transform. Used in zero padding"""
    I2 = 0
    while 2**I2 < target:
        I2 += 1
    I3 = 0
    while 3**I3 < target:
        I3 += 1
    I5 = 0
    while 5**I5 < target:
        I5 += 1
    return int(np.min([2**I2, 3**I3, 5**I5]))


def get_interpolation_factor(dt):
    """Return the automatic sinc interpolation factor: the power of two that brings the time step
    toward 0.001 s, at most 8. Interpolation captures response peaks that fall between time steps,
    which are otherwise missed at short periods."""
    return 2 ** max(min(int(np.ceil(np.log(dt / 0.001) / np.log(2))), 3), 0)


def sinc_interpolate(data, interpolation_factor):
    """Return data sinc interpolated to interpolation_factor (a power of two) times as many points,
    by zero padding its Fourier transform. Ported from Interpft in RCTC (Wang et al., MIT license,
    https://github.com/wltcwpf/RCTC)."""
    data = np.asarray(data, dtype=float)
    if interpolation_factor == 1:
        return data
    npts = len(data)
    npw2_exp = int(np.floor(np.log(npts) / np.log(2)))
    npw2 = 2**npw2_exp if 2**npw2_exp == npts else 2 ** (npw2_exp + 1)
    padded = np.concatenate([data, np.zeros(npw2 - npts)])
    m = len(padded)
    nyquist = npw2 // 2  # zero-based index of the Nyquist term
    a = np.fft.fft(padded)
    b = np.concatenate(
        [a[: nyquist + 1], np.zeros(m * interpolation_factor - m, dtype=complex), a[nyquist + 1 :]]
    )
    b[nyquist] = b[nyquist] / 2
    b[nyquist + m * interpolation_factor - m] = b[nyquist]
    y = np.real(np.fft.ifft(b) * interpolation_factor)
    return y[: npts * interpolation_factor]


def _nigam_jennings_coefficients(T, D, dt):
    """Recursion coefficients of Nigam and Jennings (1969) for an oscillator with period T and
    damping D, exact for input acceleration varying linearly over each time step."""
    omega_n = 2 * np.pi / T
    omega_d = omega_n * np.sqrt(1 - D**2)
    k = omega_n**2
    e = np.exp(-D * omega_n * dt)
    s = np.sin(omega_d * dt)
    c = np.cos(omega_d * dt)
    r = np.sqrt(1 - D**2)
    A = e * (D / r * s + c)
    B = e * (1 / omega_d * s)
    C = (
        1
        / k
        * (
            2 * D / omega_n / dt
            + e * (((1 - 2 * D**2) / omega_d / dt - D / r) * s - (1 + 2 * D / omega_n / dt) * c)
        )
    )
    Dc = 1 / k * (1 - 2 * D / omega_n / dt + e * (((2 * D**2 - 1) / omega_d / dt) * s + 2 * D / omega_n / dt * c))
    AA = -e * (omega_n / r * s)
    BB = e * (c - D / r * s)
    CC = 1 / k * (-1 / dt + e * ((omega_n / r + D / dt / r) * s + 1 / dt * c))
    DD = 1 / k / dt * (1 - e * (D / r * s + c))
    return A, B, C, Dc, AA, BB, CC, DD


def get_relative_displacement(T, D, motion, dt):
    """Return the relative displacement time series of a single-degree-of-freedom oscillator
    with period T and damping D, in the units of motion times s^2. The Nigam and Jennings (1969)
    recursion, starting at rest,

        x[j] = A x[j-1] + B v[j-1] - C a[j-1] - Dc a[j]
        v[j] = AA x[j-1] + BB v[j-1] - CC a[j-1] - DD a[j]

    is run as a second-order IIR filter after eliminating v."""
    a = np.asarray(motion, dtype=float)
    x = np.zeros(len(a))
    if len(a) < 2:
        return x
    A, B, C, Dc, AA, BB, CC, DD = _nigam_jennings_coefficients(T, D, dt)
    x[1] = -C * a[0] - Dc * a[1]
    if len(a) > 2:
        b = [-Dc, -B * DD - C + Dc * BB, -B * CC + C * BB]
        den = [1.0, -(A + BB), A * BB - B * AA]
        zi = lfiltic(b, den, y=[x[1], x[0]], x=[a[1], a[0]])
        x[2:], _ = lfilter(b, den, a[2:], zi=zi)
    return x


def _resolve_interpolation_factor(interpolation_factor, dt):
    if interpolation_factor == "auto":
        return get_interpolation_factor(dt)
    if interpolation_factor < 1 or int(interpolation_factor) & (int(interpolation_factor) - 1):
        raise ValueError("interpolation_factor must be 'auto' or a power of 2")
    return int(interpolation_factor)


def get_pseudo_acceleration_histories(T, D, motion, dt, interpolation_factor="auto"):
    """Return pseudo-acceleration time series, omega_n^2 * u(t), of single-degree-of-freedom oscillators

    Inputs:
        T = numpy array of spectral periods in seconds
        D = damping for single-degree-of-freedom oscillator
        motion = Numpy array containing an input acceleration record
        dt = time step in seconds
        interpolation_factor = 'auto' or a power of 2. The motion is sinc interpolated by this factor before
            the oscillator response is computed

    Returns a len(T) x (len(motion) * interpolation_factor) Numpy array, with time step dt / interpolation_factor.
    T = 0 returns the interpolated motion itself.
    """
    factor = _resolve_interpolation_factor(interpolation_factor, dt)
    fine = sinc_interpolate(motion, factor)
    return np.asarray([_pseudo_acceleration(Tk, D, fine, dt / factor) for Tk in np.atleast_1d(T)])


def _pseudo_acceleration(T, D, motion, dt):
    if T <= 0:
        return np.asarray(motion, dtype=float)
    return (2 * np.pi / T) ** 2 * get_relative_displacement(T, D, motion, dt)


def get_ARS(T, D, motions, dt, zeropad=None, interpolation_factor="auto"):
    """Return pseudo-acceleration response spectrum for a Numpy array containing acceleration records

    Inputs:
        T = numpy array of spectral periods in seconds
        D = damping for single-degree-of-freedom oscillator
        motions = Numpy array containing input acceleration records
        dt = time step in seconds
        zeropad = ignored, retained for backward compatibility
        interpolation_factor = 'auto' or a power of 2 for sinc interpolation of the motions

    Returns an M x len(T) Numpy array of pseudo-spectral accelerations in the units of motions.
    T = 0 returns PGA of the motions as recorded.
    """
    motions = np.atleast_2d(motions)
    T = np.atleast_1d(T)
    factor = _resolve_interpolation_factor(interpolation_factor, dt)
    Sa = np.empty((len(motions), len(T)))
    for i, motion in enumerate(motions):
        fine = sinc_interpolate(motion, factor)
        for k, Tk in enumerate(T):
            if Tk <= 0:
                Sa[i, k] = np.max(np.abs(motion))
            else:
                Sa[i, k] = np.max(np.abs(_pseudo_acceleration(Tk, D, fine, dt / factor)))
    return Sa


def _angle_bound(x, y, norms):
    """A lower bound on the peak at each of the 180 one-degree angles: the
    smallest over those angles of the peak of the 256 points farthest from
    the origin, rotated. A point closer to the origin than it can't be the
    peak at any of those angles. Zero for 512 points or fewer, which aren't
    worth pruning."""
    if len(x) <= 512:
        return 0.0
    far = np.argpartition(norms, -256)[-256:]
    theta = np.radians(np.arange(0, 180))
    return np.min(
        np.max(
            np.abs(np.outer(np.cos(theta), x[far]) + np.outer(np.sin(theta), y[far])),
            axis=1,
        )
    )


def _hull_cones(data1, data2):
    """The vertices of the convex hull of the points (data1, data2) and their
    reflections through the origin, as distance r and direction phi, and the
    range of directions [a, b] over which each vertex sets the peak: there the
    peak of the two series rotated to angle theta is r cos(theta - phi).

    Points that can't be a vertex are dropped first. _angle_bound() bounds
    the peak only at one-degree angles, and the least peak can lie between
    them, far below it for a strongly polarized orbit. The hull of the points
    kept gives a bound over all directions, its least peak, which is no more
    than the whole hull's: the points dropped that are farther from the
    origin than that are added to the hull's vertices, and the hull found
    again from those."""
    norms = np.hypot(data1, data2)
    bound = _angle_bound(data1, data2, norms)
    keep = norms >= bound
    vertices = _hull_vertices(np.column_stack([data1[keep], data2[keep]]))
    cones = _vertex_cones(vertices)
    least = _least_peak(*cones)
    if least < bound:
        missed = (norms >= least) & ~keep
        if np.any(missed):
            cones = _vertex_cones(
                _hull_vertices(
                    np.vstack(
                        [vertices, np.column_stack([data1[missed], data2[missed]])]
                    )
                )
            )
    return cones


def _least_peak(r, phi, a, b):
    """The least peak over all directions: at the end of each vertex's range
    farther from its direction."""
    return float(np.min(r * np.cos(np.maximum(np.abs(a - phi), np.abs(b - phi)))))


def _hull_vertices(points):
    """The vertices, counterclockwise, of the convex hull of points (rows
    (x, y)) and their reflections through the origin."""
    both = np.vstack([points, -points])
    try:
        # In 2-D, Qhull lists the vertices counterclockwise.
        return both[ConvexHull(both).vertices]
    except (QhullError, ValueError):
        # On a line through the origin (or all zero): the farthest point and
        # its reflection, each setting the peak over half the circle.
        far = points[np.argmax(np.hypot(points[:, 0], points[:, 1]))]
        return np.array([far, -far])


def _vertex_cones(vertices):
    """_hull_cones() of a hull's vertices, counterclockwise (_hull_vertices())."""
    r = np.hypot(vertices[:, 0], vertices[:, 1])
    phi = np.arctan2(vertices[:, 1], vertices[:, 0])
    if len(vertices) == 2:
        a, b = phi - np.pi / 2, phi + np.pi / 2
    else:
        edges = np.roll(vertices, -1, axis=0) - vertices
        # Outward normal of each counterclockwise edge
        normal = np.arctan2(-edges[:, 0], edges[:, 1])
        a, b = np.roll(normal, 1), normal
    # Each range runs counterclockwise from a to b. A vertex's own direction
    # needn't lie in its range, but every direction in the range is within
    # 90 degrees of it, so phi is taken as the copy nearest the range's middle.
    b = a + np.mod(b - a, 2 * np.pi)
    phi = phi + 2 * np.pi * np.round((0.5 * (a + b) - phi) / (2 * np.pi))
    return r, phi, a, b


def _fraction_below(levels, r, phi, a, b):
    """For each level m, the fraction of directions in which the peak is at
    most m: in each vertex's range, those where |theta - phi| >= arccos(m / r)."""
    beta = np.arccos(np.clip(np.asarray(levels, dtype=float)[:, None] / r, -1.0, 1.0))
    above = np.clip(np.minimum(b, phi + beta) - np.maximum(a, phi - beta), 0, None)
    return np.sum((b - a) - above, axis=1) / (2 * np.pi)


def _cone_percentiles(r, phi, a, b, percentiles):
    """Percentiles over all directions of the peak, given the hull's vertices
    and their ranges (_hull_cones()), found by bisection."""
    percentiles = np.asarray(percentiles, dtype=float)
    # Most at the farthest vertex, whose own direction lies in its range.
    least = _least_peak(r, phi, a, b)
    most = float(np.max(r))
    out = np.where(percentiles <= 0, least, most).astype(float)
    inner = (percentiles > 0) & (percentiles < 100)
    if np.any(inner):
        lo = np.full(np.sum(inner), least)
        hi = np.full(np.sum(inner), most)
        target = percentiles[inner] / 100
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            below = _fraction_below(mid, r, phi, a, b) < target
            lo = np.where(below, mid, lo)
            hi = np.where(below, hi, mid)
            if np.all(hi - lo <= 4 * np.finfo(float).eps * most):
                break
        out[inner] = 0.5 * (lo + hi)
    return out


def rotd_percentiles(data1, data2, percentiles):
    """Return percentiles over all angles from 0 to 180 degrees of the peak absolute value of
    data1*cos(theta) + data2*sin(theta), i.e., RotDxx for xx in percentiles

    The peak at any angle is the extent in that direction of the convex hull of the orbit
    (data1, data2) and its reflection through the origin. Between changes of hull vertex the
    peak is a cosine of the angle, so the fraction of angles below any level is computed
    exactly, and RotDxx is the level at which that fraction equals xx%, found by bisection.
    RotD0 (the nearest hull edge) and RotD100 (the farthest vertex) are exact outright. The
    result is exact over the continuous range of angles rather than over a discrete set.
    """
    data1 = np.asarray(data1, dtype=float)
    data2 = np.asarray(data2, dtype=float)
    percentiles = np.atleast_1d(percentiles)
    if len(data1) == 0 or not np.any(data1) and not np.any(data2):
        return np.zeros(len(percentiles))
    return _cone_percentiles(*_hull_cones(data1, data2), percentiles)


def get_rotD_ARS(T, D, motion1, motion2, dt, zeropad=None, RotD=50, N_rot_angle=None, interpolation_factor="auto"):
    """Return RotD pseudo-acceleration response spectrum for two Numpy arrays containing orthogonal horizontal acceleration records

    Inputs:
        T = numpy array of spectral periods in seconds
        D = damping for single-degree-of-freedom oscillator
        motion1 = Numpy array containing input acceleration record for motion 1
        motion2 = Numpy array containing input acceleration record for motion 2
        dt = time step in seconds
        zeropad = ignored, retained for backward compatibility
        RotD = Numpy array containing desired RotD values (e.g., RotD = np.asarray([10,50,90]))
        N_rot_angle = ignored, retained for backward compatibility. RotD values are computed
            exactly over all angles using the convex hull method (see rotd_percentiles)
        interpolation_factor = 'auto' or a power of 2 for sinc interpolation of the motions

    Returns a len(T) x len(RotD) Numpy array of pseudo-spectral accelerations in the units of the motions.
    T = 0 returns RotD of PGA of the motions as recorded.
    """
    if len(motion2) != len(motion1):
        raise ValueError("Both motions must have same length")
    T = np.atleast_1d(T)
    RotD = np.atleast_1d(np.asarray(RotD, dtype=float))
    factor = _resolve_interpolation_factor(interpolation_factor, dt)
    fine1 = sinc_interpolate(motion1, factor)
    fine2 = sinc_interpolate(motion2, factor)
    Sa = np.empty((len(T), len(RotD)))
    for k, Tk in enumerate(T):
        if Tk <= 0:
            Sa[k] = rotd_percentiles(motion1, motion2, RotD)
        else:
            Sa[k] = rotd_percentiles(
                _pseudo_acceleration(Tk, D, fine1, dt / factor),
                _pseudo_acceleration(Tk, D, fine2, dt / factor),
                RotD,
            )
    return Sa


def get_response_spectrum(**kwargs):
    """Function that reads keyword arguments and calls get_ARS or get_rotD_ARS function

    Keyword arugments:
        motions = numpy array containing acceleration time series (required)
        dt = time step for acceleration records (required)
        T = numpy array of spectral periods (default = NGAWest2 spectral periods)
        D = damping for single-degree-of-freedom oscillator
        RotD = RotD percentile or list of percentiles (e.g., 50 or [0, 50, 100]). If specified,
            motions must contain two orthogonal horizontal components
        interpolation_factor = 'auto' (default) or a power of 2 for sinc interpolation of the motions
        zeropad = ignored, retained for backward compatibility
        N_rot_angle = ignored, retained for backward compatibility

    Returns pseudo-spectral accelerations in the units of motions: an M x len(T) array, or a
    len(T) x len(RotD) array if RotD is specified.
    """
    if "motions" not in kwargs:
        print("you must specify motions")
        return
    if "dt" not in kwargs:
        print("you must specify dt")
        return
    T = np.atleast_1d(np.asarray(kwargs.get("T", get_ngawest2_T()), dtype=float))
    motions = np.asarray(kwargs["motions"], dtype=float)
    if len(motions.shape) == 1:
        motions = motions.reshape((1, len(motions)))
    D = kwargs.get("D", 0.05)
    dt = kwargs["dt"]
    interpolation_factor = kwargs.get("interpolation_factor", "auto")
    if "RotD" in kwargs:
        if len(motions) != 2:
            raise ValueError("RotD requires motions to contain two orthogonal horizontal components")
        return get_rotD_ARS(
            T, D, motions[0], motions[1], dt, RotD=kwargs["RotD"], interpolation_factor=interpolation_factor
        )
    else:
        return get_ARS(T, D, motions, dt, interpolation_factor=interpolation_factor)
