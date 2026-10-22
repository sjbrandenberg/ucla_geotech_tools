# random_field

This function generates the random field using both Cholesky Decomposition method and Kriging. Details are explained by Yang et al. (2022).  

Yang, Y., Wang, P., and Brandenberg, S. J. (2022) "An algorithm for generating spatially correlated random fields using Cholesky decomposition and ordinary kriging." Computers and Geotechnics, 147. [[official version](https://doi.org/10.1016/j.compgeo.2022.104783)] [[open access](https://authors.elsevier.com/a/1f13c,63b%7EyYS2)]

The nodes are divided into Cholesky nodes and kriging nodes. Each node is visited in order, and becomes a kriging node if it is within Dhchol of a Cholesky node, and otherwise a Cholesky node. The field is generated at the Cholesky nodes by Cholesky decomposition of their correlation matrix, and at the kriging nodes by ordinary kriging from the Cholesky nodes. Kriging estimates are smoother than the random field, so an independent normal residual with the kriging variance is added at each kriging node, which gives the field the variance sigma^2 at every node.

Kriging reproduces a smooth (Gaussian) correlation structure from relatively widely spaced Cholesky nodes, but a rough (exponential) correlation structure requires closely spaced Cholesky nodes. The default values of Dhchol (0.15 theta_z for the Gaussian model and 0.07 theta_z for the exponential model) limit the error in the correlation between pairs of nodes to about 0.03. Larger values of Dhchol are faster but less accurate.

## Installation
```bash
pip install ucla_geotech_tools
```

## Functions
```bash
get_random_field(**kwargs): main function for generating the random field
get_chol_krig_nodes(x, z, Dhchol): divide nodes into Cholesky nodes and kriging nodes, and return their indices and coordinates
get_C(x, z, theta, model): correlation matrix among nodes
get_C_chol_krig(xchol, zchol, xkrig, zkrig, theta, model): correlation matrix between Cholesky nodes and kriging nodes
get_chol_field(C, mu, sigma): random field at Cholesky nodes
get_krig_field(C, C_chol_krig, chol_field): ordinary kriging estimates at kriging nodes
```

## Input Parameters for main function
```bash
x: x coordinates (list or array of x coordinates)
z: z coordinates (list or array of z coordinates)
mu: mean of random field
sigma: standard deviation of random field
theta_z: vertical scale of fluctuation
model: (optional) 0 for Gaussian correlation (default), 1 for exponential (Markovian) correlation
theta_x: (optional) horizontal scale of fluctuation (default theta_z)
Dhchol: (optional) maximum distance from a kriging node to the nearest Cholesky node, after scaling x by theta_z/theta_x (default 0.15*theta_z for Gaussian, 0.07*theta_z for exponential)
seed: (optional) seed or numpy random Generator for a reproducible field. If not specified, numpy's global random generator is used, so np.random.seed() also gives reproducible fields.
```

## Output Data
```bash
output : an array of random field values at the nodes, in the same order as x and z.
```

## Usage
```python
import ucla_geotech_tools.random_field as ugtrf
import numpy as np

x = np.linspace(0.0,40,641)
z = np.linspace(0.0,10,161)
x,z = np.meshgrid(x,z)
x = x.reshape(np.prod(x.shape),) # create a vector of X coordinates
z = z.reshape(np.prod(z.shape),) # create a vector of Y coordinates

RF = ugtrf.get_random_field(x=x,z=z,mu=50,sigma=10,theta_z=1.25,seed=1)
```
