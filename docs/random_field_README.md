# random_field

This function generates the random field using both Cholesky Decomposition method and Kriging

## Installation
```bash
pip install ucla_geotech_tools.random_field
```

## Functions
```bash
sort_nodes(x,z): sort the nodes to better separate cholesky nodes from kriging nodes; return sorted node coordinates x and z
get_random_field(x,z,theta_z,mu,sigma): main function for generating the random field
get_chol_krig_nodes(x,z,Dhchol): group cholesky nodes and kriging nodes separately and return the index of cholesky and kriging nodes
get_COV(xvchol,zvchol,theta_z,model,sigma_chol): get covariance matrix
get_C_chol_krig(xvchol, zvchol, xvkrig, zvkrig, theta_z, model, sigma_chol, sigma_krig): get covariance matrix between cholesky nodes and kriging nodes
get_chol_field(COV, mu_chol): generate the random field for cholskey nodes
get_krig_field(COV, C_chol_krig, chol_field): generate the random field for kriging nodes
```

## Input Parameters for main function
```bash
x: x coordinates (list or array of x coordinates)
z: z coordinates (list or array of z coordinates)
mu: mean of random field
sigma: standard deviation of random field
theta_z: vertical scale of fluctuation
model: (optional) default is 0 for Guassian correlation matrix. Other intergers for Markovian correlation matrix
```

## Output Data
```bash
output : an array of values, and it has the same of order of the sorted nodes.
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
x,z = ugtrf.sort_nodes(x,z)      # sort nodes

RF = ugtrf.get_random_field(x=x,z=z,mu=50,sigma=10,theta_z = 1.25)
```
