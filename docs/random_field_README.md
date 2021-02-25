# random_field

This function generates the random field using both Cholesky Decomposition method and Kriging

## Installation
```bash
pip install ucla_geotech_tools.random_field
```

```bash
Note: Please make sure you have the package "pyKrige" installed before installing this package.
You can try "pip install pykrige" if you haven't done so!
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
import ucla_geotech_tools.random_field as rf
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0.0,40,41)
y = np.linspace(0.0,40,41)
X,Y = np.meshgrid(x,y)
X = X.reshape(np.prod(X.shape),1) # create a vector of X coordinates
Y = Y.reshape(np.prod(Y.shape),1) # create a vector of Y coordinates
kwargs = {"X" : X, "Y" : Y, "corr_func" : 'Gaussian',"beta_h" : 8,"beta_v" : 8,"mu" : 40,"c1" : 100,"M" : 2,"beta":0.6} 
data = rf.random_field(**kwargs) # output data will have a dimension of (len(X),M)
data1 = data[:,0].reshape(41,41) # only use the first column of values to plot the random field
plt.imshow(data1)
plt.colorbar()
plt.show()
```

## Example
The Jupyter notebook example illustrates how to use the random_field package. Two examples are given in this Jupyter notebook, where random fields are generated using Gaussian and Markovian coorelation function, respectively. Feel free to modify the example to suit your own problem. Right click the link and use "Save link as".

[random_field.ipynb](https://raw.githubusercontent.com/sjbrandenberg/ucla_geotech_tools/main/random_field/random_field.ipynb)
