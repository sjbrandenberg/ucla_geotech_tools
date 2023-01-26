# auto_fchp

auto_fchp selects a high-pass corner frequency for a single component ground motion acceleration record, acc. The algorithm fits a polynomial to the filtered displacement record, and iterates on the high-pass corner frequency until the ratio of the amplitude of the polynomial fit to that of the displacement record is equal to a specified target. The algorithm follows these steps:

1.	Subtract the mean from acc
2.	Apply a Tukey window 
3.	Compute the Fourier transform, <img src="https://render.githubusercontent.com/render/math?math=F_{acc}">, and frequency vector, f
4.	Select a trial high-pass corner frequency, <img src="https://render.githubusercontent.com/render/math?math=\hat{f}_{chp}">
5.	Filter the record using an acausal Butterworth filter defined by Eq. 1, where u is an index counter over the number of frequency components
6.	Compute the Fourier coefficients of the displacement record, <img src="https://render.githubusercontent.com/render/math?math=F_{disp}">, using Eq. 2
7.	Compute the displacement time series, disp, by computing the inverse Fourier transform of <img src="https://render.githubusercontent.com/render/math?math=F_{disp}">
8.	Fit a polynomial of a desired order, <img src="https://render.githubusercontent.com/render/math?math=disp_{fit}">, to disp
9.	Compute the value of the error function, E, defined by Eq. 3 where target is the desired value of the ratio of the amplitude of the polynomial fit to that of the displacement

Equation 1:  
  
<img src="https://render.githubusercontent.com/render/math?math=filter_u = \frac{1}{\sqrt{1+\left(\frac{\hat{f}_{chp}}{f_u}\right)^{2\cdot order}}}">

Equation 2:  
  
<img src="https://render.githubusercontent.com/render/math?math=Fdisp_u = \frac{Facc_u \cdot filter_u}{-\left(2\pi f_u\right)^2}"> 

Equation 3:  
  
<img src="https://render.githubusercontent.com/render/math?math=E = \frac{\left|disp_{fit}\right|}{\left|disp\right|} - target">  

## Installation  
```python
pip install ucla_geotech_tools.auto_fchp
```

## Function
```python
get_fchp(**kwargs)
maxabs(double[:] vx)
filtered_Facc(complex[:] Facc, double[:] freq, double fc, double order)
get_vel(double[:] freq, complex[:] Facc)
get_disp(double[:] freq, complex[:] Facc)

Note: the [:] operator indicates that these functions represent arrays using typed memoryview objects for efficiency.
https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html
```

## Input parameters
```
### get_fchp(**kwargs) -> double fchp

    Compute the high-pass corner frequency that renders a tolerable displacement time series after high-pass filtering an
    acceleration record. Tolerable is defined based on the ratio of the computed displacement trace to that of a polynomial
    that is fit to the displacement trace.
    
    Note: **kwargs is an arbitrary keyword argument that indicates many different types of input parameters may be input
    
    Required inputs
    double dt -- time step in seconds
    numpy.ndarray acc -- acceleration timeseries formatted as a Numpy NDarray object
    
    Optional inputs
    target -- desired ratio of polynomial fit to displacement amplitude (default 0.02).
    tol -- tolerance for ratio of polynomial fit to displacement amplitude (default 0.001).
    poly_order -- order of polynomial fit to displacement record (default 6)
    maxiter -- maximum number of iterations (default 30)
    fchp_min -- minimum permissible value of fchp in Hz (default 0.001)
    fchp_max -- maximum permissible value of fchp in Hz (default 0.5)
    filter_order -- high pass Butterworth filter order (default 5)
    tukey_alpha -- percent of time series to apply Tukey window at beginning and end of record (default 0.20)
    filter_type -- type of filter to apply (default 1)
                   filter_type = 1 is a frequency-domain Butterworth filter
                   filter_type = 0 is a finite impulse response filter designed using Obspy  (https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.filter.html)
    
    example command: fchp = get_fchp(dt=dt, acc=acc)

### maxabs(double[:] vx) -> double output

    Compute the maximum of the absolute value of real valued array vx.
    vx = array of double-precision values.
    Note, the [:] notation indicates that vx is handled as a typed 
    memoryview (https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html).

### filtered_Facc(complex[:] Facc, double[:] freq, double fc, double order) -> complex[:] output
    Compute the Fourier coefficients after applying a high-pass Butterworth filter
    complex[:] Facc = array of complex-valued Fourier coefficients obtained from running a fast Fourier transform.
    double[:] freq = array of frequency values associated with the Fourier transform
    double fc = high pass corner frequency for Butterworth filter
    double order = order for Butterworth filter
    Note, the [:] notation indicates that an array is handled as a typed memoryview
    memoryview (https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html).
    Returns a complex-valued typed memoryview object containing filtered Fourier coefficients.
    The typed memoryview can be converted to a Numpy array like this: output_array = np.asarray(output)

### get_vel(double[:] freq, complex[:] Facc) -> double[:] vel 
    Compute the velocity time series using frequency-domain integration of the acceleration trace Fourier coefficients
    double[:] freq = array of frequency values associated with the Fourier transform
    complex[:] Facc = array of complex-valued Fourier coefficients obtained from running a fast Fourier transform.
    Note, the [:] notation indicates that an array is handled as a typed memoryview
    memoryview (https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html).
    Returns a typed memoryview object containing the velocity time series.
    The typed memoryview can be converted to a Numpy array like this: vel_array = np.asarray(vel)

### get_disp(double[:] freq, complex[:] Facc) -> double[:] disp
    Compute the displacement time series using frequency-domain integration of the acceleration trace Fourier coefficients
    double[:] freq = array of frequency values associated with the Fourier transform
    complex[:] Facc = array of complex-valued Fourier coefficients obtained from running a fast Fourier transform.
    Note, the [:] notation indicates that an array is handled as a typed memoryview
    memoryview (https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html).
    Returns a typed memoryview object containing the displacement time series.
    The typed memoryview can be converted to a Numpy array like this: vel_array = np.asarray(vel)
```

## Subroutines

The get_fchp function calls a number of subroutines that are defined using cdef in Cython, and are therefore not accessible via the pip installable package. However, you may download and modify the source code if you would like to utilize these functions directly. These functions include maxabs(vx), which returns the maximum of the absolute value of a vector vx, filters_Facc(Facc, freq, fc, order), whic returns Fourier coefficients for the filtered version of Facc, get_vel(freq, Facc), which returns the Fourier coefficients for velocity given a frequency vector (freq) and Fourier coefficients for the acceleration (Facc), get_disp(freq, Facc), which returns Fourier coefficents for displacement, and get_residual(time, disp, target, poly_order), which returns the residual defined by Eq. 3.
