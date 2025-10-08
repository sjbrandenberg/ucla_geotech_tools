"""random_field is a tool for efficiently creating realizations of two-dimensional random fields by combining Cholesky decomposition of a covariance matrix with kriging at intermediate nodes.

Details of the method are described by Yang and Brandenberg (2022).

Yang, Y., Wang, P., and Brandenberg, S. J. (2022) "An algorithm for generating spatially correlated random fields using Cholesky decomposition and ordinary kriging." Computers and Geotechnics, 147

"""
import numpy as np

def get_C(x, z, theta, model):
    """Returns a correlation matrix for an input array containing x and z coordinates, a correlation length, and the type of correlation model (0 = Gaussian, 1 = exponential)
    
    Args:
        x: array containing x coordinates
        z: array containing z coordinates
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    dist = np.sqrt((x - x[:, np.newaxis]) ** 2.0 + (z - z[:, np.newaxis]) ** 2.0)
    if(model == 0):
        C = np.exp(-np.pi * (dist / theta) ** 2.0)
    else:
        C = np.exp(-2.0 * np.abs(dist) / theta)
    
    return C

def get_C_chol_krig(xchol, zchol, xkrig, zkrig, theta, model):
    """Returns correlation matrix between the cholesky decomposition nodes and the kriging nodes.

    Args:
        xchol: array containing x coordinates of cholesky nodes
        zchol: array containing z coordinates of cholesky nodes
        xkrig: array containing x coordinates of kriging nodes
        zkrig: array containing z coordinates of kriging nodes
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    dist = np.sqrt((xchol - xkrig[:, np.newaxis]) ** 2.0 + (zchol - zkrig[:, np.newaxis]) ** 2.0)
    if(model == 0):
        C = np.exp(-np.pi * (dist / theta) ** 2.0)
    else:
        C = np.exp(-2.0 * np.abs(dist) / theta)

    return C

def get_chol_krig_nodes(x, z, dist_max):
    """Divides nodes into a set for cholesky decomposition and a set for kriging.
    
    Args:
    x: x coordinates of all nodes
    z: z coordinates of all nodes
    dist_max: maximum distance by which cholesky nodes may be separated
    """
    ind = np.arange(len(x))
    ind_chol = np.asarray([0], dtype=int)
    ind_krig = np.asarray([], dtype=int)

    while(len(ind_chol) + len(ind_krig) < len(ind)):
        mask = ~np.isin(ind, np.hstack((ind_chol, ind_krig)))
        dist = np.sqrt((x[ind_chol[-1]] - x) ** 2 + (z[ind_chol[-1]] - z) ** 2)
    
        # If we have found an isolated set where all points are beyond max_dist, pick a new point that is not already selected
        if(np.min(dist[mask]) >= dist_max):
            ind_chol = np.append(ind_chol, ind[mask][0])
        else:
            ind_chol = np.append(ind_chol, np.argwhere(dist == np.max(dist[mask][dist[mask] < dist_max]))[0])
        ind_krig = np.append(ind_krig, ind[(dist < dist_max) & ~np.isin(ind, np.hstack((ind_chol, ind_krig)))])
        
    return (ind_chol, x[ind_chol], z[ind_chol], ind_krig, x[ind_krig], z[ind_krig])

def get_chol_field(C, mu, sigma):
    """Function that uses cholesky decomposition to return a random realization at the cholesky nodes

    Args:
        C: correlation matrix at cholesky nodes
        mu: mean of random field
        sigma: standard deviation of random field    
    """
    Nchol = len(C)
    rv = np.random.normal(loc=0, scale=1, size=Nchol)
    L_chol = np.linalg.cholesky(C+0.0001*np.identity(Nchol))
    chol_field = mu + sigma*L_chol@rv

    return chol_field

def get_krig_field(C, C_chol_krig, chol_field):
    """Returns values of random field at kriging nodes

    Args:
        C: correlation matrix at cholesky nodes
        C_chol_krig: correlation matrix between cholesky nodes and kriging nodes
        chol_field: realization of random field at cholesky nodes
    """
    Nchol = len(C)
    Nkrig = C_chol_krig.shape[0]
    C_chol_krig = np.vstack((C_chol_krig.T, np.ones(Nkrig)))
    C = np.hstack((C,np.transpose([np.ones(Nchol)])))
    C = np.vstack((C,np.ones(Nchol+1)))
    C[Nchol,Nchol] = 0.0
    krig = np.linalg.solve(C, C_chol_krig)
    W = krig[0:Nchol, 0:Nkrig]
    krig_field = np.transpose(W)@chol_field
    return krig_field

def get_random_field(**kwargs):
    """Function that returns composite random field from cholesky decomposition combined with kriging.
    
    Keyword Args:
        x: x coordinates of all nodes (required)
        z: z coordinates of all nodes (required)
        theta_z: correlation length in z-direction (required)
        mu: mean of random field (required)
        sigma: standard deviation of random field (required)
        model: correlation model, 0 = Gaussian 1 = exponential (default 0)
        theta_x: correlation length in x-direction (default theta_z)
        Dhchol: maximum length by which cholesky nodes are allowed to be separated (default np.sqrt(2.0*np.pi)*theta_z/4.0)
    """
    required = ['x', 'z', 'theta_z', 'mu', 'sigma']
    for r in required:
        if r not in kwargs.keys():
            print(r + ' must be included like this "' + r + '= value"')
            return
    x = np.asarray(kwargs['x'])
    z = np.asarray(kwargs['z'])
    theta_z = kwargs['theta_z']
    mu = kwargs['mu']
    sigma = kwargs['sigma']
    model = kwargs.get('model',0)
    theta_x = kwargs.get('theta_x',theta_z)
    Dhchol = kwargs.get('Dhchol',np.sqrt(2.0*np.pi)*theta_z/4.0)
       
    x = x*theta_z/theta_x # perform coordinate transformation
    chol_ind, xvchol, zvchol, krig_ind, xvkrig, zvkrig  = get_chol_krig_nodes(x, z, Dhchol) # get cholesky and kriging nodes
    C = get_C(xvchol, zvchol, theta_z, model) # get correlation matrix
    chol_field = get_chol_field(C, mu, sigma) # get cholesky decomposition field
    C_chol_krig = get_C_chol_krig(xvchol, zvchol, xvkrig, zvkrig, theta_z, model) # get covariance matrix between kriging and cholesky points
    krig_field = get_krig_field(C, C_chol_krig, chol_field) # get kriging field
    rf = np.zeros(len(x)) # initialize random field array
    for i, ind in enumerate(chol_ind):
        rf[ind] = chol_field[i] # populate random field at cholesky nodes
    for i, ind in enumerate(krig_ind):
        rf[ind] = krig_field[i] # populate random field at kriging nodes

    return rf