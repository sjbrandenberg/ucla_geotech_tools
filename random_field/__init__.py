def random_field(X,Y,corr_func,beta_h,beta_v,mu,c1,M,beta = 0.6):
    """
    This function generates the random field using both Cholesky Decomposition method and Kriging
    X : X coordinates
    Y : Y coordinates
    corr_func : Correlation function, "Gaussian" or "Markovian"
    beta_h : horizontal scale of fluctuation
    beta_v : vertical scale of fluctuation
    mu : mean of random field
    c1 : variance of random field
    M : number of iterations preferred
    beta : a scalar defining sample coarseness
    """
    import numpy as np
    from scipy.spatial import distance_matrix
    import pykrige.kriging_tools as kt
    from pykrige.ok import OrdinaryKriging
    
    data = np.zeros((len(X),M))
    Z = np.arange(0,len(X)).reshape(len(X),1)
    coords1= np.concatenate((Z,X,Y),axis=1)
    coords2= np.concatenate((X,Y),axis=1)
    coords11 = coords1[np.lexsort(coords2.T)]
    Z = coords11[:,0].astype(int)
    coords = coords2[np.lexsort(coords2.T)]
    X = coords[:,0]
    Y = coords[:,1]
    X_chol = []
    Y_chol = []
    Z_chol = []
    X_chol.append(X[0])
    Y_chol.append(Y[0])
    Z_chol.append(Z[0])
    for x,y,z in zip(X,Y,Z):
        min_dist = 10**30
        for x_chol,y_chol in zip(X_chol,Y_chol):
            if ((x-x_chol)**2+(y-y_chol)**2)**0.5 < min_dist:
                min_dist = ((x-x_chol)**2+(y-y_chol)**2)**0.5
        if min_dist >= beta*beta_h:
            X_chol.append(x)
            Y_chol.append(y)
            Z_chol.append(z)
            
    if beta == 0.0:
        Z_chol = np.delete(Z_chol,0,0)
        X_chol = np.delete(X_chol,0)
        Y_chol = np.delete(Y_chol,0)
    X_chol =np.reshape(X_chol,(len(X_chol),1))
    Y_chol =np.reshape(Y_chol,(len(Y_chol),1))
    coords_chol = np.concatenate((X_chol*beta_v/beta_h,Y_chol),axis=1)
    if beta != 0.0 and len(Z_chol) != len(Z):
        Z_krig = np.asarray([ele for ele in Z if ele not in Z_chol])
        X_krig =X[Z_krig].reshape(len(X[Z_krig]),1)
        Y_krig = Y[Z_krig].reshape(len(Y[Z_krig]),1)
        coords_krig = np.concatenate((X_krig*beta_v/beta_h,Y_krig),axis=1)
    sigma_ln = np.sqrt(np.log(1+c1/(mu**2)))
    mu_ln = np.log(mu)-1/2*sigma_ln**2

    dist1 = distance_matrix(coords_chol,coords_chol)
    if corr_func == 'Gaussian':
        C1 = sigma_ln**2*np.exp(-np.pi*(dist1/beta_v)**2)
    if corr_func == 'Markovian':
        C1 = sigma_ln**2*np.exp(-2*(dist1/beta_v))
        
    L_chol = np.linalg.cholesky(C1+0.0001*np.identity(len(coords_chol)))

    if beta == 0.0:
        for i in range(M):
            e = np.random.normal(0,1,len(coords_chol))
            data_ln = mu_ln + L_chol @ e
            data1 = np.exp(data_ln)
            data[Z_chol,i] = data1
    if beta != 0.0:
        if len(Z_chol) == len(Z):
            for i in range(M):
                e = np.random.normal(0,1,len(coords_chol))
                data_ln = mu_ln + L_chol @ e
                data1 = np.exp(data_ln)
                data[Z_chol,i] = data1
        if len(Z_chol) != len(Z):
            for i in range(M):
                e = np.random.normal(0,1,len(coords_chol))
                data_ln = mu_ln + L_chol @ e
                data1 = np.exp(data_ln)
                if corr_func == 'Gaussian':
                    OK = OrdinaryKriging(coords_chol[:,0],coords_chol[:,1],data_ln,variogram_model='gaussian',variogram_parameters=[sigma_ln**2,beta_v,0])
                if corr_func == 'Markovian':
                    OK = OrdinaryKriging(coords_chol[:,0],coords_chol[:,1],data_ln,variogram_model='exponential',variogram_parameters=[sigma_ln**2,beta_v,0])
                    
                z, ss = OK.execute('points', coords_krig[:,0],coords_krig[:,1])

                data[Z_chol,i] = data1
                data[Z_krig,i] = np.exp(z)
    return data