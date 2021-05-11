## ipyconsol

ipyconsol is a nonlinear implicit finite difference solver for one-dimensional consolidation of compressible soil with secondary compression. ipyconsol is the Python version of the Javascript package iConsol.js published by Brandenberg (2017).

### Installation  
```bash
pip install ucla_geotech_tools.ipyconsol
```

### Input Parameters

#### General input parameters.
```bash
qo = initial vertical effective stress at the top of the layer  
tol = convergence tolerance (optional, default = 1.0e-8)
gammaw = unit weight of water (optional, default = 9.81)
pa = atmospheric pressure (optional, default = 101.325)
drainagetype = 0 double drainage, 1 single drainage through the top, or 2 single drainage through the bottom (optional, default = 0)
```
  
#### Soil properties (uniform soil layer). All inputs are scalars
```bash
H = layer thickness
N = number of elements
Cc = virgin compression index
Cr = recompression index
sigvref = vertical effective stress of reference point on normal consolidation line
esigvref = void ratio of reference point on normal consolidation line
Gs = specific gravity of solids
Ck = coefficient of permeability variation
kref = hydraulic conductivity of reference point on e-logk curve
ekref = void ratio of reference point on e-logk curve
Ca = secondary compression index
tref = reference time corresponding to normal consolidation line
dsigv = total stress increment applied to top of soil layer
ocrvoidratiotype = 0 constant OCR, 1 constant eo, 2 constant maximum past pressure
ocrvoidratio = value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype
ru = initial excess pore pressure, prior to addition of dsigv (optional, default = 0.0)
```
#### Soil properties (layered profile). All inputs are lists or numpy arrays
```bash
depth = node locations
Cc = virgin compression index values
Cr = recompression index values
sigvref = vertical effective stresses for reference point on normal consolidation line
esigvref = void ratios for reference point on normal consolidation line
Gs = specific gravities of solids
Ck = coefficients of permeability variation
kref = hydraulic conductivities for reference point on e-logk curve
ekref = void ratios for reference point on e-logk curve
Ca = secondary compression indices
tref = reference times corresponding to normal consolidation line
dsigv = total stress increments applied at each depth
ocrvoidratiotype = 0 constant OCR, 1 constant eo, 2 constant maximum past pressure
ocrvoidratio = value of OCR, eo, or maximum past pressure, depending on value of ocrvoidratiotype
ru = initial excess pore pressures, prior to addition of dsigv (optional, default = 0.0)
```
#### Time input parameters (constant load increment)
```bash
Ntime = number of time increments
tmax = maximum time value
```
#### Time input parameters (time-dependent loading sequence)
```bash
time = list or array of time values
loadfactor = corresponding load factor values. Stress increment is load factor multiplied by dsigv.
```

### Usage

```python
import ucla_geotech_tools.ipyconsol as pcl

pcl.compute(*args, **kwargs) # returns dictionary containing depth values, 'z', pore pressures, 'u', vertical effective stresses 'sigv', and void ratios 'e' 

pcl.get_inputs(*args, **kwargs) # reads inputs and returns [depth, Cc, Cr, sigvref, esigvref, Gs, kref, ekref, Ck, Ca, tref, qo, dsigv, ocrvoidratiotype, 
                                # ocrvoidratio, ru, time, loadfactor, gammaw, tol, pa, drainagetype] 

pcl.get_initial(depth, Cc, Cr, sigvref, esigvref, Gs, kref, ekref, Ck, Ca, tref, qo, dsigv, ocrvoidratiotype, ocrvoidratio, ru, time, loadfactor, gammaw, tol, pa, drainagetype)
# returns dictionary containing initial void ratio 'eo', vertical effective stress 'sigvo', vertical effective stress with ru 'sigvo', and final effective stress 'sigvf'

pcl.Logspace(tmin, tmax, Ntime, tstart) # returns time vector evenly distributed in log space between tmin and tmax

pcl.get_ktest(etest,ekref,kref,Ck) # returns hydraulic conductivity for specified void ratio, etest

pcl.get_avtest(elast, sigvlast, sigvtest, eref, sigvref, Cc, Cr) # returns coefficient of compressibility

pcl.get_etest(etest, elast, eref, alpha, tref, sigvtest, sigvlast, sigvref, Cc, avtest, utest, ulast, dt)
# returns trial void ratio

pcl.get_residual(ulast, utest, elast, etest, klast, ktest, avtest, zlast, ztest, sigvlast, sigvtest, gammaw, Ca, Cc, sigvref, esigvref, double dt, tref, drainagetype, N)
# returns array of residual values for trial pore pressure solution

pcl.get_utest2(N, klast,  ktest,  zlast,  ztest,  avtest,  elast,  etest,  ulast,  utest,  sigvlast,  sigvtest,  Res, dt,  Cc,  Cr,  Ca, gammaw,  ekref,  kref,  Ck,  tref,  esigvref,  sigvref,  dsigv, drainagetype, pa)
# returns trial pore pressure solution using Newton Raphson iteration

pcl.get_utest(klast, ktest, avtest, zlast, ztest, ulast, sigvlast, sigvtest, elast, etest, Ca, Cc, sigvref, esigvref, dt, tref, drainagetype, N, gammaw)
# returns initial pore pressure solution using linear theory

```

### Example

The Jupyter notebook example illustrates how to use the ipyconsol package. Four different cases are analyzed: 1. Uniform soil, constant load, 2. Layered soil, constant load, 3. Uniform soil, time-dependent load, and 4. Layered soil, time-dependent load. Feel free to modify the example to suit your own problem. Right click the link and use "Save link as".

[ipyconsol.ipynb](https://github.com/sjbrandenberg/ucla_geotech_tools/raw/main/ipyconsol/ipyconsol.ipynb)

### Web version
A version of this software is available [here](https://www.uclageo.com/Consolidation) for uniform soil, constant load and [here](https://www.uclageo.com/Consolidaiton2) for layered soil, time-varying load.

### References
Brandenberg, S.J. (2017). "iConsol.js: A javascript implicit finite difference code for nonlinear consolidation and secondary compression." International Journal of Geomechanics, 17(6). [[journal](http://ascelibrary.org/doi/abs/10.1061/%28ASCE%29GM.1943-5622.0000843)] [[eScholarship](https://escholarship.org/uc/item/0wh3q8jh)]
