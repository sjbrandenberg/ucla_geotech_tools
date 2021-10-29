# random_field

This function generates the random field using both Cholesky Decomposition method and Kriging

## Installation
```bash
pip install ucla_geotech_tools.random_field
```

## Functions
```bash
sort_nodes: sort the nodes to better separate cholesky nodes from kriging nodes
get_random_field: main function for generating the random field
get_C: get correlation matrix
get_chol_krig_nodes: group cholesky nodes and kriging nodes separately
get_C_chol_krig: get correlation matrix between cholesky nodes and kriging nodes
get_chol_field: generate the random field for cholskey nodes
get_krig_field: generate the random field for kriging nodes
```

## Input Parameters for main function
```bash
x: x coordinates (list or array of x coordinates)
z: z coordinates (list or array of z coordinates)
mu: mean of random field
sigma: standard deviation of random field
theta_v: vertical scale of fluctuation
model: (optional) default is 0 for Guassian correlation matrix. Other intergers for Markovian correlation matrix
```

## Output Data
```bash
output : an array of values, and it has the same of order of the sorted nodes.
```

## Usage
```python
import ucla_geotech_tools.random_field as ugtrf
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0.0,40,641)
z = np.linspace(0.0,10,161)
x,z = np.meshgrid(x,z)
x = x.reshape(np.prod(x.shape),) # create a vector of X coordinates
z = z.reshape(np.prod(z.shape),) # create a vector of Y coordinates

x,z = ugtrf.sort_nodes(x,z)      # sort nodes

RF = ugtrf.get_random_field(x=x,z=z,mu=50,sigma=10,theta_z = 1.25)

f,ax = plt.subplots(figsize=(10, 5))
cntr = ax.tricontourf(x,z,RF,20)
f.colorbar(cntr, ax=ax,shrink=0.5,aspect = 5)
ax.set_aspect('equal')
ax.set_frame_on(False)
```

## Example
The Jupyter notebook example illustrates how to use the random_field package. Two examples are given in this Jupyter notebook, where random fields are generated using Gaussian and Markovian coorelation function, respectively. Feel free to modify the example to suit your own problem. Right click the link and use "Save link as".

[random_field.ipynb](https://github.com/sjbrandenberg/ucla_geotech_tools/blob/main/random_field/random_field.ipynb)
