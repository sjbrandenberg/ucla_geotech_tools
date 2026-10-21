# response_spectrum

response_spectrum computes a pseudo-acceleration response spectrum from an input ground motion or set of ground motions sampled at a constant frequency. It is also capable of computing rotated response spectra for two horizontal components (e.g., RotD50). Calculations for computing a response spectrum are performed in the time domain by:

1. sinc interpolating the input motion(s) to a finer time step, by zero padding the Fourier transform. By default the interpolation factor is the power of two that brings the time step toward 0.001 s, at most 8. Without interpolation, peaks of the oscillator response that fall between time steps are missed, and short-period spectral accelerations are too low.
2. computing the relative displacement {math}`u` of a single-degree-of-freedom (SDOF) oscillator with natural period {math}`T` and damping {math}`D` using the recursion of Nigam and Jennings (1969), which is exact for input acceleration varying linearly over each (interpolated) time step. The recursion is evaluated as an IIR filter using scipy.signal.lfilter.
3. computing the pseudo-spectral acceleration {math}`PSA = \omega_n^2 |u|_{max}`
4. repeating steps 2 and 3 for a range of user-specified periods. {math}`T = 0` returns PGA of the motion as recorded.

The interpolation and recursion are ported from RCTC (Wang et al., https://github.com/wltcwpf/RCTC, MIT license). Results match the NGA ground motion database processing code to machine precision.

A time-domain solution was selected over a frequency-domain solution because the frequency domain approach requires zero padding long enough for the free vibration of the longest-period oscillator to decay (approximately {math}`\ln(1000)/(D \omega_n)`, or 440 s for {math}`T` = 20 s and {math}`D` = 0.05) to prevent the response from wrapping around to the beginning of the record. With sufficient padding and interpolation the two approaches agree within 0.01%, but the time-domain approach is approximately five times faster.

RotD values are computed exactly over all rotation angles from 0 to {math}`\pi` using a convex hull method, rather than over a discrete set of angles:

1. Compute the pseudo-acceleration time series {math}`\omega_n^2 u_1(t)` and {math}`\omega_n^2 u_2(t)` of an SDOF oscillator with a specific T and D for the two orthogonal horizontal motions, {math}`motion1` and {math}`motion2`.
2. The peak of the rotated series {math}`\omega_n^2 [u_1 \cos\theta + u_2 \sin\theta]` at any angle {math}`\theta` is the extent in that direction of the convex hull of the orbit {math}`(u_1, u_2)` and its reflection through the origin. Only the hull vertices can set a peak, so the remaining points are discarded.
3. Between changes of hull vertex, the peak is a cosine function of {math}`\theta`, so the fraction of angles at which the peak is below any level can be computed exactly. RotDxx is the level at which that fraction equals xx%, found by bisection. RotD0 (the nearest hull edge) and RotD100 (the farthest hull vertex) are exact.

Computing RotD over a discrete set of angles (e.g., 1 degree increments) can underestimate the minimum, so RotD0 of a strongly polarized response can be substantially in error.

## Signal conditioning

response_spectrum does not perform any signal conditioning, which is assumed to be performed prior to invoking the response_spectrum command. I recommend performing the following steps prior to computing a response spectrum:

1. Assign a Tukey window with {math}`\alpha` = 0.05, where {math}`\alpha` denotes the percentage of the signal at the beginning and end of the record for the window ramp up/down
2. Subtract the weighted mean from the signal, where the weights are equal to the Tukey window coefficients
3. High-pass filter. You can try the ucla_geotech_tools.auto_fchp package for this.

## Installation  
```python
pip install ucla_geotech_tools
```

## Functions
```python
get_response_spectrum(**kwargs)
get_ngawest2_T()
rotd_percentiles(data1, data2, percentiles)
```

### get_response_spectrum(**kwargs)

**Keyword Args:**  

| parameter | type | description | required | default |
|-----------|------|-------------|----------|---------|   
| motions  | numpy array, dtype = float | M x N numpy array containing acceleration data, where M is the number of motions and N is the number of data points in each motion. All motions must have N data points. | x | |
| dt | float | time step | x | |
| D | float | damping | | 0.05 |
| T | numpy array, dtype = float | array of natural periods | | array used for NGAwest2 project |
| interpolation_factor | 'auto' or int | sinc interpolation factor, a power of 2. 'auto' brings the time step toward 0.001 s, at most a factor of 8 | | 'auto' |
| RotD | float or list of floats | RotD percentile(s), e.g., 50 or [0, 50, 100]. If specified, motions must be a 2 x N array containing two orthogonal horizontal components | | |
| zeropad | int | ignored; retained for backward compatibility | | |
| N_rot_angle | int | ignored; retained for backward compatibility | | |

**Returns:** pseudo-spectral accelerations in the units of motions, as an M x len(T) numpy array, or a len(T) x len(RotD) numpy array if RotD is specified.

### get_ngawest2_T

Returns numpy array containing natural periods utilized for NGAWest2 response spectra

### rotd_percentiles(data1, data2, percentiles)

Returns the percentiles over all angles from 0 to {math}`\pi` of the peak absolute value of data1*cos(theta) + data2*sin(theta), computed exactly with the convex hull method described above. Useful for computing RotD values of any pair of orthogonal time series (e.g., PGA or PGV).

### Example Commands

The example below creates 20 harmonic motions with various frequencies and computes their response spectra. Additional examples can be found in the tests folder of the GitHub repository.

```python
N = 2**14
dt = 0.005
Nmotion = 20
time = np.linspace(0,dt*N,N)
f = np.linspace(0.5,8.0,Nmotion)
motions = np.sin(2*np.pi*f[:,np.newaxis]*time)
window = sp.signal.windows.tukey(N,0.20)
mean_motions = np.average(motions,weights=window,axis=1)
motions = window*(motions - mean_motions[:,np.newaxis])
T = rs.get_ngawest2_T()
Sa = rs.get_response_spectrum(motions=motions,D=0.05,T=T,dt=dt)
