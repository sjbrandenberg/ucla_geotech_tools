# random_field

This function generates the random field using both Cholesky Decomposition method and Kriging

## Installation
```bash
pip install ucla_geotech_tools.random_field
```

## Input Parameters
```bash
X : X coordinates (list or array of X coordinates)
Y : Y coordinates (list or array of Y coordinates)
corr_func : Correlation function, "Gaussian" or "Markovian"
beta_h : horizontal scale of fluctuation
beta_v : vertical scale of fluctuation
mu : mean of random field
c1 : variance of random field
M : number of iterations preferred
beta : a scalar defining sample coarseness (default value is 0.6)
```

## Output Data
```bash
output : an array of values, whose Dimension is (len(X),M). Every column is one realization.
```

## Usage
```python
%matplotlib inline
from ucla_geotech_tools.random_field import random_field
import matplotlib.pyplot as plt
import numpy as np

x1 = np.linspace(0.0,40,41)
y1 = np.linspace(0.0,40,41)
X,Y = np.meshgrid(x1,y1)
X = X.reshape(np.prod(X.shape),1) # create a vector of X coordinates
Y = Y.reshape(np.prod(Y.shape),1) # create a vector of Y coordinates
kwargs = {"X" : X, "Y" : Y, "corr_func" : 'Gaussian',"beta_h" : 8,"beta_v" : 8,"mu" : 40,"c1" : 100,"M" : 2,"beta":0.6} 
data = random_field(**kwargs) # output data will have a dimension of (len(X),M)
data1 = data[:,0].reshape(41,41) 
plt.imshow(data1)
plt.colorbar()
plt.show()
```

## Change Log
```bash
0.0.1(02/02/2021)
-First Release
```
