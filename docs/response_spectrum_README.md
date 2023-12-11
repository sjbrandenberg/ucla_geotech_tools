# response_spectrum

response_spectrum computes a psuedo-acceleration response spectrum from an input ground motion or set of ground motions sampled at a constant frequency. It is also capable of computing rotated response spectra for two horizontal components (e.g., RotD50). Calculations for computing a response spectrum are performed in the frequency domain by:

1. taking the Fourier transform of the input motion(s)
2. convolving with the transfer function for the resonse of a single-degree-of-freedom (SDOF) oscillator with natural period $T$ and damping $$
4. computing the peak acceleration amplitude for the SDOF oscillator
5. repeating steps 2 and 3 for a range of user-specified periods

Calculations for RotD50 are performed by:

1. Computing the response of an SDOF oscillator with a specific T and D for the two orthogonal horizontal motions, $motion1$, and $motion2$.
2. Setting _**N_rot_angle**_ different rotation angles between 0 and $\pi$
3. Computing convolved motion for each angle as motion_rot_conv = motion1_conv*cos(angle) - motion2_conv*sin(angle)
4. Compute the maximum of each value of motion_rot_conv. This is called motion_rot_conv_max, and there is one value per angle.
5. Compute the specified percentile value of the motion_rot_conv_max.

## Signal conditioning

response_spectrum does not perform any signal conditioning, which is assumed to be performed prior to invoking the response_spectrum command. I recommend performing the following steps prior to computing a response spectrum:

1. Assign a Tukey window with $\alpha$ = 0.05, where $\alpha$ denotes the percentage of the signal at the beginning and end of the record for the window ramp up/down
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
```

### get_response_spectrum(**kwargs)

**Keyword Args:**  

| parameter | type | description | required | default |
|-----------|------|-------------|----------|---------|   
| motions  | numpy array, dtype = float | M x N numpy array containing acceleration data, where M is the number of motions and N is the number of data points in each motion. All motions must have N data points. | x | |
| dt | float | time step | x | |
| D | float | damping | | 0.05 |
| T | numpy array, dtype = float | array of natural periods | | array used for NGAwest2 project |
| zeropad | int | 1: apply zero padding to speed up FFT operation, 0: do not pad with zeros | | 1 |

### get_ngawest2_T

Returns numpy array containing natural periods utilized for NGAWest2 response spectra

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
Fmotions = np.fft.rfft(motions)
freq = np.fft.rfftfreq(N,dt)
T = rs.get_ngawest2_T()
Sa = rs.get_response_spectrum(motions=motions,D=0.05,T=T,dt=dt,zeropad=0)
