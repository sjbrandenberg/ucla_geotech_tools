# cython: language_level=3

"""random_field is a tool for efficiently creating realizations of two-dimensional random fields by combining Cholesky decomposition of a covariance matrix with kriging at intermediate nodes.

Details of the method are described by Yang and Brandenberg (2022).

Yang, Y., Wang, P., and Brandenberg, S. J. (2022) "An algorithm for generating spatially correlated random fields using Cholesky decomposition and ordinary kriging." Computers and Geotechnics, 147

"""
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, sqrt, pow, pi, abs

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
def get_C(double[:] x, double[:] z, double theta, int model):
    """Returns a correlation matrix for an input array containing x and z coordinates, a correlation length, and the type of correlation model (0 = Gaussian, 1 = exponential)
    
    Args:
        x: array containing x coordinates
        z: array containing z coordinates
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    cdef double [:, :] C = np.zeros((len(x), len(z)), dtype='float64')
    cdef int i, j
    cdef double dist
    for i in range(len(x)):
        for j in range(i,len(z),1):
            dist = sqrt(pow(x[i]-x[j],2.0) + pow(z[i]-z[j],2.0))
            if(model==0):
                C[i][j] = exp(-pi*(pow(dist/theta,2.0)))
                C[j][i] = C[i][j]
            else:
                C[i][j] = exp(-2*abs(dist)/theta)
                C[j][i] = C[i][j]
    return np.asarray(C)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
def get_C_chol_krig(double[:] xchol, double[:] zchol, double[:] xkrig, double[:] zkrig, double theta, int model):
    """Returns correlation matrix between the cholesky decomposition nodes and the kriging nodes.

    Args:
        xchol: array containing x coordinates of cholesky nodes
        zchol: array containing z coordinates of cholesky nodes
        xkrig: array containing x coordinates of kriging nodes
        zkrig: array containing z coordinates of kriging nodes
        theta: correlation length
        model: 0 = Gaussian, 1 = exponential
    """
    cdef double [:, :] C = np.zeros((len(xchol), len(xkrig)), dtype='float64')
    cdef int i, j
    cdef double dist
    for i in range(len(xchol)):
        for j in range(len(xkrig)):
            dist = sqrt(pow(xchol[i]-xkrig[j],2.0) + pow(zchol[i]-zkrig[j],2.0))
            if(model==0):
                C[i][j] = exp(-pi*(pow(dist/theta,2.0)))
            else:
                C[i][j] = exp(-2*abs(dist)/theta)
    return np.asarray(C)
       
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
def get_chol_krig_nodes(double[:] x, double[:] z, double dist_max):
    """Divides nodes into a set for cholesky decomposition and a set for kriging.
    
    Args:
    x: x coordinates of all nodes
    z: z coordinates of all nodes
    dist_max: maximum distance by which cholesky nodes may be separated
    """
    
    cdef double[:] x2z2 = np.zeros(len(x))
    cdef int i, ind
    cdef double xmax, zmax, xzmax
    chol_ind = [0]
    xvchol = [x[0]]
    zvchol = [z[0]]
    krig_ind = []
    xvkrig = []
    zvkrig = []
    cdef double dist
    cdef int krig = 0
    for i in range(1, x.size):
        krig = 0
        for k in range(len(xvchol)):
            if(sqrt((x[i]-xvchol[k])*(x[i]-xvchol[k]) + (z[i]-zvchol[k])*(z[i]-zvchol[k])) < dist_max):
                krig_ind.append(i)
                xvkrig.append(x[i])
                zvkrig.append(z[i])
                krig = 1
                break
        if(krig==0):
            chol_ind.append(i)
            xvchol.append(x[i])
            zvchol.append(z[i])
    return [np.asarray(chol_ind), np.asarray(xvchol), np.asarray(zvchol), np.asarray(krig_ind), np.asarray(xvkrig), np.asarray(zvkrig)]


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
    Nkrig = C_chol_krig.shape[1]
    C_chol_krig = np.vstack((C_chol_krig, np.ones(Nkrig)))
    C = np.hstack((C,np.transpose([np.ones(Nchol)])))
    C = np.vstack((C,np.ones(Nchol+1)))
    C[Nchol,Nchol] = 0.0
    krig = np.linalg.solve(C, C_chol_krig)
    W = krig[0:Nchol, 0:Nkrig]
    krig_field = np.transpose(W)@chol_field
    return krig_field

@cython.boundscheck(False)
@cython.wraparound(False)
def sort_nodes(double[:] x, double[:] z):
    """Sorts nodes by spatial position to improve efficiency of algorithm to identify cholesky nodes.
    
    Args:
        x: x coordintes for all nodes
        z: z coordintes for all nodes
    """
    cdef int i
    values = []
    for i in range(len(x)):
        values.append((i+1, sqrt(pow(x[i]-x[0],2.0) + pow(z[i]-z[0],2.0)), x[i], z[i]))
    sort_mat = np.sort(np.array(values, dtype=[('ind', int), ('dist', float), ('xv', float), ('zv', float)]), order='dist')
    x = sort_mat['xv']
    z = sort_mat['zv']
    return([np.asarray(x), np.asarray(z)])

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
    ind_old = np.arange(1,len(x),1)
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
