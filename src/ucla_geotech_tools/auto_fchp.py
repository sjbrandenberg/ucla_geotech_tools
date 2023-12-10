import scipy as sp
import numpy as np

def get_residual1(fchp, *args):
    """Return the residual defined as disp_fit/disp_filt - target.

    Args:
        fchp = High-pass corner frequency
        Fdisp = Numpy array containing Fourier coefficients of displacement record
        time = Numpy array containing time values
        poly_order = Order of polynomial to fit to filtered displacement
        target = Desired ratio of polynomial fit amplitude to displacement signal amplitude
        freq = Numpy array containing frequency values

    Computed variables:
        filt = Coefficients of 5th order Butterworth filter
        Fdisp_filt = Numpy array containing Fourier coefficients of filtered displacement record
        disp_filt = Numpy array containing filtered displacement record
        coef = Coefficients of polynomial fit to displacement record
        disp_fit = Numpy array containing polynomial fit ordinates at each time
    """
    Fdisp = args[0]
    time = args[1]
    poly_order = args[2]
    target = args[3]
    freq = args[4]
    filter_order = args[5]
    filt = np.zeros(len(Fdisp))
    filt[freq > 0] = 1 / np.sqrt(1 + (fchp / freq[freq > 0]) ** (2 * filter_order))
    Fdisp_filt = Fdisp * filt
    disp_filt = np.fft.irfft(Fdisp_filt)
    coef = np.polyfit(time[0 : len(disp_filt)], disp_filt, poly_order)
    disp_fit = np.polyval(coef, time)
    return np.max(np.abs(disp_fit)) / np.max(np.abs(disp_filt)) - target


def get_residual2(fchp, *args):
    """Return the residual defined as disp_init/disp - target.

    Args:
        fchp = High-pass corner frequency
        Fdisp = Numpy array containing Fourier coefficients of displacement record
        time = Numpy array containing time values
        disp_time = Duration of initial portion of record before p-wave arrival
        disp_target = Desired ratio of amplitude of initial portion of displacement record to amplitude of displacement record
        freq = Numpy array containing frequency values

    Computed variables:
        filt = Coefficients of 5th order Butterworth filter
        Fdisp_filt = Numpy array containing Fourier coefficients of filtered displacement record
        disp_filt = Numpy array containing filtered displacement record
    """
    Fdisp = args[0]
    time = args[1]
    disp_time = args[2]
    disp_target = args[3]
    freq = args[4]
    filter_order = args[5]
    filt = np.zeros(len(Fdisp))
    filt[freq > 0] = 1 / np.sqrt(1 + (fchp / freq[freq > 0]) ** (2 * filter_order))
    Fdisp_filt = Fdisp * filt
    disp_filt = np.fft.irfft(Fdisp_filt)
    return (
        np.max(np.abs(disp_filt[time <= disp_time])) / np.max(np.abs(disp_filt))
        - disp_target
    )


def get_fchp(**kwargs):
    """Return the high-pass corner frequency value that stabilizes the displacement of a ground motion record.

    The high-pass corner frequency, fchp, is selected to satisfy two objectives:
        1. A polynomial is fit to the filtered displacement record, and the amplitude of the polynomial fit
        divided by the amplitude of the filtered displacement record is within a tolerance of a specified target.
        2. The amplitude of the initial portion of the displacement record before the p-wave arrival divided by
        the amplitude of the filtered displacement record is within a tolerance of a specified target.

    The value of fchp is computed numerically using the Scipy package optimize.ridder(). Two initial values
    of fchp are provided that bracket the root, and the Scipy package solves for the root that minimizes the
    objective function. The tolerance is on the value of fchp (i.e., xtol in optimize.ridder)

    A value of fchp is obtained first to satisfy objective 1. Objective 2 is optional, and if it is applied,
    the value of fchp from objective 1 is checked against objective 2. If the residual from objective 2 is
    less than or equal to zero, fchp is set equal to the value from objective 1. Otherwise, iterations are
    performed to satisfy objective 2.

    Keyword Args:
        dt (float): time step in seconds (required)
        acc (numpy array, dtype=float): acceleration array in L/s/s (required)
        target: desired ratio of polynomial fit amplitude to displacement amplitude (default=0.02)
        tol: tolerance for fchp (default=0.001)
        poly_order: order of polynomial to fit to filtered displacement record (default=6)
        maxiter: maximum number of iterations (default=30)
        fchp_min: minimum frequency in Hz to consider for fchp (default=0.001)
        fchp_max: maximum frequency in Hz to consider for fchp (default=0.5)
        tukey_alpha: alpha parameter for scipy.signal.windows.tukey (default=0.05)
        apply_disp_ratio: flag to indicate whether to apply objective 2 (default=0)
        disp_ratio_time: duration in seconds of the initial portion of the record before the p-wave arrival (default=30)
        disp_ratio_target: target ratio of amplitude of initial portion of displacement record to amplitude of displacement record
    """
    options = [
        "dt",
        "acc",
        "target",
        "tol",
        "poly_order",
        "maxiter",
        "fchp_min",
        "fchp_max",
        "filter_order",
        "tukey_alpha",
        "apply_disp_ratio",
        "disp_ratio_time",
        "disp_ratio_target",
    ]

    for key, value in kwargs.items():
        if key not in options:
            print(
                key
                + " is not a valid argument. Please see documentation. Using default values for all parameters that are not specified."
            )
            return
    if('dt' not in kwargs):
        print('You must specify dt')
        return
    if('acc' not in kwargs):
        print('You must specify acc')
        return
    dt = kwargs['dt']
    acc = kwargs['acc']
    target = kwargs.get("target", 0.02)
    tol = kwargs.get("tol", 0.001)
    poly_order = kwargs.get("poly_order", 6)
    maxiter = kwargs.get("maxiter", 30)
    fchp_min = kwargs.get("fchp_min", 0.001)
    fchp_max = kwargs.get("fchp_max", 0.5)
    filter_order = kwargs.get("filter_order", 5.0)
    tukey_alpha = kwargs.get("tukey_alpha", 0.05)
    apply_disp_ratio = kwargs.get("apply_disp_ratio", 0)
    disp_ratio_time = kwargs.get("disp_ratio_time", 30.0)
    disp_ratio_target = kwargs.get("disp_ratio_target", 0.05)

    # Signal pre-processing. Subtract weighted mean, apply tukey window, define time, take FFT, compute Fdisp
    if len(acc) % 2 != 0:
        acc = acc[0:-1]
    window = sp.signal.windows.tukey(len(acc), alpha=tukey_alpha)
    meanacc = np.average(acc, weights=window)
    acc = window * (acc - meanacc)
    time = np.linspace(0, dt * len(acc), len(acc))
    Facc = np.fft.rfft(acc)
    freq = np.fft.rfftfreq(len(acc), dt)
    Fdisp = np.zeros(len(Facc), dtype="complex")
    Fdisp[freq > 0] = Facc[freq > 0] / -((2 * np.pi * freq[freq > 0]) ** 2)

    # Solve fchp for objective 1
    R0 = get_residual1(fchp_min, Fdisp, time, poly_order, target, freq, filter_order)
    R2 = get_residual1(fchp_max, Fdisp, time, poly_order, target, freq, filter_order)
    if np.sign(R0) != np.sign(R2):
        fchp = sp.optimize.ridder(
            get_residual1,
            fchp_min,
            fchp_max,
            args=(Fdisp, time, poly_order, target, freq, filter_order),
            xtol=tol,
            maxiter=maxiter,
        )
    elif np.sign(R0) < 0:
        fchp = fchp_min
    elif np.sign(R2) > 0:
        fchp = fchp_max

    # Solve fchp for objective 2 if apply_disp_ratio is not zero
    if apply_disp_ratio == 0:
        return fchp
    R0 = get_residual2(
        fchp, Fdisp, time, disp_ratio_time, disp_ratio_target, freq, filter_order
    )
    R2 = get_residual2(
        fchp_max, Fdisp, time, disp_ratio_time, disp_ratio_target, freq, filter_order
    )
    if np.sign(R0) != np.sign(R2):
        fchp = sp.optimize.ridder(
            get_residual2,
            fchp,
            fchp_max,
            args=(Fdisp, time, disp_ratio_time, disp_ratio_target, freq, filter_order),
            xtol=tol,
            maxiter=maxiter,
        )
    elif np.sign(R0) < 0:
        fchp = fchp
    elif np.sign(R2) > 0:
        fchp = fchp_max

    return fchp
