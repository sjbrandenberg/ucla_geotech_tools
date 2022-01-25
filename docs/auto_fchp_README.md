# auto_fchp

auto_fchp selects a high-pass corner frequency for a single component ground motion acceleration record, acc. The algorithm fits a polynomial to the filtered displacement record, and iterates on the high-pass corner frequency until the ratio of the amplitude of the polynomial fit to that of the displacement record is equal to a specified target. The algorithm follows these steps:

1.	Subtract the mean from acc
2.	Apply a Tukey window 
3.	Compute the Fourier transform, <img src="https://render.githubusercontent.com/render/math?math=F_{acc}">, and frequency vector, f
4.	Select a trial high-pass corner frequency, <img src="https://render.githubusercontent.com/render/math?math=\hat{f}_{chp}">
5.	Filter the record using an acausal Butterworth filter defined by Eq. 1
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
```

## Input parameters
```
Required keys:
    
    acc = A numpy array or Python list defining acceleration time series.
    dt = time step
        
Optional keys:
    
    target = target ratio of polynomial fit amplitude to displacement amplitude (default = 0.02)
    tol = tolerance for fchp (default = 0.001 Hz)
    filter_order = Butterworth filter order (default = 5)
    poly_order = order of polynomial fit (default = 6)
    fchp_min = minimum value of range of permissible corner frequencies (default = 0.001 Hz)
    fchp_max = maximum value of range of permissible corner frequencies (default = 0.5 Hz)
    maxiter = maximum number of iterations (default = 30)
    tukey_alpha = alpha parameter defining width of Tukey window (default = 0.2)
        
Note: **kwargs defines a Python dictionary in the form "key=value". The key must be defined exactly as shown above, 
and the value is the data quantity associated with the key. For example:

acc = [0.00, 0.01, -0.002, ..., 0.04, -0.001]
dt = 0.005
get_response_spectra(acc=acc, dt=dt)

In this case, default values will be assigned to all but the two required parameters.
```

## Subroutines

The get_fchp function calls a number of subroutines that are defined using cdef in Cython, and are therefore not accessible via the pip installable package. However, you may download and modify the source code if you would like to utilize these functions directly. These functions include maxabs(vx), which returns the maximum of the absolute value of a vector vx, filters_Facc(Facc, freq, fc, order), whic returns Fourier coefficients for the filtered version of Facc, get_vel(freq, Facc), which returns the Fourier coefficients for velocity given a frequency vector (freq) and Fourier coefficients for the acceleration (Facc), get_disp(freq, Facc), which returns Fourier coefficents for displacement, and get_residual(time, disp, target, poly_order), which returns the residual defined by Eq. 3.
