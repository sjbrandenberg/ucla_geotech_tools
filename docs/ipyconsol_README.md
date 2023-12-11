# ipyconsol

ipyconsol is a nonlinear implicit finite difference solver for one-dimensional consolidation of compressible soil with secondary compression. 
ipyconsol is the Python version of the Javascript package iConsol.js published by Brandenberg (2017).
This code is an extension of Brandenberg (2017) because it the soil can be non-uniform and the loading can vary in time.
For spatially variable soil, material constants are entered as arrays, where the array length is equal to the number of nodes.
For time-varying loading, a loadfactor array is specified and has a length equal to the time array.

Brandenberg, S. J. (2017) "iConsol. js: JavaScript implicit finite-difference code for nonlinear consolidation and secondary compression." International Journal of Geomechanics, 17(6)

## Installation  
```bash
pip install ucla_geotech_tools
```
## Requirements
numpy >= 1.22  
C compiler (e.g., gcc, Visual Studio Build Tools)

## Functions  

```python
compute(**kwargs)
get_initial(depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype)
get_inputs(**kwargs)
```

## Variables

There are two options for soil inputs (uniform and non-uniform) and two options for loading (constant and time-varying).

soil property options:
    option 1: soil is uniform (i.e., material parameters are constants)
    option 2: soil is non-uniform (i.e., material parameters are arrays, with array lengths equal to number of nodes)

load stage options:
    option A: constant load applied at time=0
    option B: time-varying load

The options may be combined in any manner (e.g., [option 1, option A], [option 1, option B], [option 2, option A], [option 2, option B])
Uniform soil has constant model coefficients, and either void ratio, OCR, or maximum past pressure is constant.

### Keyword Args:

| parameter | type | description | required | default |
|-----------|------|-------------|----------|---------|
| ```N``` | int | number of elements | required for option 1 | |
| ```H``` | float | thickness of soil layer | required for option 1 | |
| ```depth``` | numpy array, dtype = float | array of depth values | required for option 2 | |
| ```Ntime``` | float | number of time steps | required for option A | |
| ```tmax``` | float | maximum time value | required for option A | |
| ```time``` | numpy array, dtype = float | array of time values | required for option B | |
| ```loadfactor``` | numpy array, dtype = float | array of load factor values where stress increment = loadfactor*dsigv | required for option B | |
| ```Cc``` | float or numpy array | virgin compression index | required | | 
| ```Cr``` | float or numpy array | recompression index | required | | 
| ```sigvref``` | float or numpy array | vertical effective stress for a point on the normal consolidation line | required | | 
| ```esigvref``` | float or numpy array | void ratio for point on normal consolidation line corresponding to esigvref | required | | 
| ```Gs```| float or numpy array | specific gravity of solids | required | |
| ```kref``` | float or numpy array | hydraulic conductivity for a point on the e-logk line | required | |
| ```ekref``` | float or numpy array | void ratio for point on the e-logk like corresponding to kref | |
| ```Ck``` | float or numpy array | coefficient of permeability variation defined as slope of e-logk relationship de/dlogk | required | |
| ```Ca``` | float or numpy array | coefficient of secondary compression | required | |
| ```tref``` | float or numpy array | time associated with normal consolidation line following Bjerrum's time-line concept | required | |
| ```qo``` | float | initial vertical effective stress at top of soil | required | |
| ```dsigv``` | float or numpy array | change in total stress | required | |
| ```ocrvoidratiotype``` | float or numpy array | 0 OCR, 1 eo, 2 maximum past pressure | required | |
| ```ocrvoidratio``` | float or numpy array | value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype | required | |
| ```gammaw``` | float | unit weight of water | | 9.81 |
| ```tol``` | float | convergence tolerance | | 1.0e-8 |
| ```pa``` | float | atmospheric pressure | | 101.325 |
| ```drainagetype``` | float | 0 = double-drained, 1 = single-drained through top, 2 = single-drained through bottom | | 0 |

### Example Commands

The example command below is for a uniform soil layer and constant loading. For other examples, see the test_ipyconsol.py file in the GitHub repository.

```python
import ucla_geotech_tools.ipyconsol as pcl

output = pcl.compute(
    Cc=5.729,
    Cr=0.532,
    sigvref=50.0,
    esigvref=8.763,
    Gs=1.686,
    kref=1.0e-7,
    ekref=8.344,
    Ck=200.0,
    Ca=0.0,
    tref=3070.0,
    qo=10.0,
    dsigv=400.0,
    ocrvoidratiotype=0,
    ocrvoidratio=1.5,
    drainagetype=0,
    ru=0.0,
    pa=101.325,
    tol=1.0e-8,
    N=100,
    H=4.35,
    Ntime=500,
    tmax=35532000.0
)
```
