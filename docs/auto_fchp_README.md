# auto_fchp

auto_fchp selects a high-pass corner frequency $f_{chp}$ for a single component ground motion acceleration record, $acc$. The algorithm uses two criteria to select $f_{chp}$:

**criterion 1:** The amplitude of a polynomial fit to the displacement record must be a specified multiple of the amplitude of the displacement record  

**criterion 2:** The amplitude of the initial portion of the displacement record before the p-wave arrival must be less than or equal to a specified multiple of the amplitude of the displacement record

The algorithm first selects $f_{chp}$ to satisfy criterion 1. It then checks criterion 2, and if necessary, increases $f_{chp}$ to satify criterion 2. Criterion 2 is optional.

## criterion 1
The algorithm fits a polynomial of user-specified order to the filtered displacement record, and iterates on $f_{chp}$  until the ratio of the amplitude of the polynomial to that of the displacement record is equal to a specified target. The acceleration record is pre-processed following these steps:

1.	Define a Tukey window
2.	Subtract the weighted mean from $acc$, where the weights are equal to the Tukey window coefficients 
3.	Compute the Fourier transform, $F_{acc}$, and frequency array, $f$
4.  Compute the Fourier coefficients of the displacement record using Eq. 1

After pre-processing, a check is performed to determine whether the optimal value of $f_{chp}$ lies between $f_{chp,min}$ and $f_{chp,max}$. The check follows these steps:

5. Set $f_{chp}$ equal to $f_{chp,min}$
6. Compute the filtered Fourier displacement coefficients $Fdisp_{filt}$ using Eq. 2.
7. Compute the filtered displacement record by taking the inverse Fourier transform of $Fdisp_{filt}$
8. Compute the residual $R1$ using Eq. 3
9. Repeat steps 5 through 8 using $f_{chp,max}$ instead of $f_{chp,min}$
10. If the sign of the residuals are equal, the root is not bracketed. Return $f_{chp,max}$ if the sign is positive, and return $f_{chp,min}$ if the sign is negative.

Finally, if the root is bracketed, use Ridders' method to find the value of $f_{chp}$ that satisfies Eq. 3. In this case, the Scipy package scipy.optimize.ridder is utilized.

11. If the root is bracketed, solve for $f_{chp}$ using scipy.optimize.ridder()

Equation 1:  
  
$Fdisp = \frac{Facc}{-\left(2\pi f\right)^2}$

Equation 2:  
  
$Fdisp_{filt} = \frac{Fdisp}{\sqrt{1+\left(\frac{f_{chp}}{f_u}\right)^{2\cdot order}}}$

Equation 3:  
  
$R1 = \frac{\left|disp_{fit}\right|}{\left|disp\right|} - target$

## criterion 2 (optional)
If the user chooses to use criterion 2, the algorithm checks the ratio of the initial portion of the displacement record with duration of _**disp_ratio_time**_ and compares it with the amplitude of the filtered displacement record. If the ratio is larger than _**disp_ratio_target**_, the algorithm iterates on $f_{chp}$ until the ratio is equal to _**disp_ratio_target**_. Criterion 2 follows these steps:

12. Extract the initial portion of the displacmeent record $disp_{init}$
13. Using Eq. 4, compute the criterion 2 residual $R2$ for the value of $f_{chp}$ obtained for criterion 1.
14. If $R2 \leq tol$, return $f_{chp}$
15. If $R2 > tol$, follow the logic from criterion 1 using the scipy.optimize.ridder package to solve for $f_{chp}$ 

Equation 4:  
  
$R2 = \frac{\left|disp_{init}\right|}{\left|disp\right|} - target$

## Installation  
```python
pip install ucla_geotech_tools
```

## Requirements
numpy >= 1.22  
scipy >= 1.8

## Function
```python
get_fchp(**kwargs)
get_residual1(fchp *args)
get_residual2(fchp *args)
```

## Input parameters
### get_fchp(**kwargs)
Return the high-pass corner frequency value that stabilizes the displacement of a ground motion record.  
  
**Keyword Args:**  

| parameter | type | description | required | default |
|-----------|------|-------------|----------|---------|
|```dt```   | float | time step in seconds (required)  |  x  |  |
|```acc```  |numpy array, dtype=float | acceleration array in L/s/s (required) |  x  |  |
|```target```| float | desired ratio of polynomial fit amplitude to displacement amplitude | | 0.02 |
|```tol```| float | tolerance for fchp | | 0.001 |  
|```poly_order```| int | order of polynomial to fit to filtered displacement record | | 6 |  
|```maxiter```| int | maximum number of iterations | | 30 |  
|```fchp_min```| float | minimum frequency in Hz to consider for fchp| | 0.001 |  
|```fchp_max```| float | maximum frequency in Hz to consider for fchp | | 0.5 |  
|```tukey_alpha```| float | alpha parameter for scipy.signal.windows.tukey | | 0.05 |  
|```apply_disp_ratio```| int | flag to indicate whether to apply criterion 2 | | 0 |  
|```disp_ratio_time```| float | duration in seconds of the initial portion of the record before the p-wave arrival | | 30.0 |  
|```disp_ratio_target```| float | target ratio of amplitude of initial portion of displacement record to amplitude of displacement record | | 0.05 |
  
**Example Commands:**  
```python
fchp = get_fchp(dt=dt, acc=acc)

fchp = get_fchp(dt=dt,acc=acc,target=0.02,tol=0.001,poly_order=6,maxiter=30,fchp_min=0.001,fchp_max=0.5,filter_order=5.0,tukey_alpha=0.05,apply_disp_ratio=1,disp_ratio_time=2,disp_ratio_target=0.02)
```

### get_residual1(fchp, *args):
Return the residual defined as disp_fit/disp_filt - target.  
  
** Args:**  

| parameter | type | description |
|-----------|------|-------------|
|```fchp```| float | High-pass corner frequency |
|```Fdisp```| numpy array, dtype = complex | Fourier coefficients of displacement record |
|```time```| numpy array, dtype = float | time array in seconds |
|```poly_order```| int |Order of polynomial to fit to filtered displacement |
|```target```| float | desired ratio of polynomial fit amplitude to displacement signal amplitude |
|```freq```| numpy array, dtype = float| frequency array |

### get_residual2(fchp, *args):
Return the residual defined as disp_init/disp - target.

**Args:**  

| parameter | type | description |
|-----------|------|-------------|
|```fchp```| float | High-pass corner frequency |
|```Fdisp```| numpy array, dtype = complex | Fourier coefficients of displacement record |
|```time``` | numpy array, dtype = float | time array in seconds |
|```disp_time```| float | duration of initial portion of record before p-wave arrival |
|```disp_target```| float | desired ratio of amplitude of initial portion of displacement record to amplitude of displacement record |
|```freq```| numpy array, dtype = float | frequency array |
