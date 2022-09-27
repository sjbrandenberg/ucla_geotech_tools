from libc.math cimport sin, cos, pi, atan, pow, log, sqrt, abs, exp
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef int next_fast_len(target):
    cdef int I2=0
    while(pow(2.0,I2)<target):
        I2=I2+1
    cdef int I3=0
    while(pow(3.0,I3)<target):
        I3=I3+1
    cdef int I5=0
    while(pow(5.0,I5)<target):
        I5=I5+1
    return(np.min([2.0**I2,3.0**I3,5.0**I5]))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef get_ARS(double[:] T, double D, double[:,:] motions, double dt, int zeropad):
    cdef int LEN = len(motions[0])
    for motion in motions:
        if(len(motion) !=LEN):
            return('Length of all motions must be the same')
    cdef int NFFT
    if(zeropad==1):
        NFFT = next_fast_len(LEN)
    else:
        NFFT = LEN
    cdef double[:] freq = np.fft.rfftfreq(NFFT,dt)
    cdef int i, j, u
    cdef double complex[:] Fmotion = np.empty(len(freq),dtype='complex')
    cdef double complex[:] FConvolvedMotion = np.empty(len(freq), dtype='complex')
    cdef double[:] ConvolvedMotion = np.empty(NFFT,dtype='float64')
    cdef double[:,:] Sa = np.empty([len(motions),len(T)],dtype='float64')
    cdef double Saval
    for i in range(len(motions)):
        Fmotion = np.fft.rfft(np.asarray(motions[i]),NFFT)
        for j in range(len(T)):
            for u in range(len(freq)):
                FConvolvedMotion[u] = Fmotion[u]/(1.0+2.0j*D*freq[u]*T[j]-freq[u]*freq[u]*T[j]*T[j])
            ConvolvedMotion = np.fft.irfft(FConvolvedMotion)
            Saval = 0.0
            for u in range(len(ConvolvedMotion)):
                if(abs(ConvolvedMotion[u])>Saval):
                    Sa[i][j] = abs(ConvolvedMotion[u])
                    Saval = Sa[i][j]
    return(np.asarray(Sa))
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def get_rotD_ARS(double[:] T, double D, double[:] motion1, double[:] motion2, double dt, int zeropad, double [:] RotD, int N_rot_angle):
    cdef int LEN = len(motion1)
    if(len(motion2) != LEN):
        return('Both motions must have same length')
    cdef int NFFT
    if(zeropad==1):
        NFFT = next_fast_len(LEN)
    else:
        NFFT = LEN
    cdef double [:] theta = np.linspace(0,np.pi,N_rot_angle)
    cdef double [:,:] rot = np.asarray([np.cos(theta), -np.sin(theta)])
    cdef double[:] freq = np.fft.rfftfreq(NFFT,dt)
    cdef int i, k, u
    cdef double complex[:] Fmotion1 = np.empty(len(freq),dtype='complex')
    cdef double complex[:] Fmotion2 = np.empty(len(freq),dtype='complex')
    cdef double complex[:] FConvolvedMotion1 = np.empty(len(freq), dtype='complex')
    cdef double complex[:] FConvolvedMotion2 = np.empty(len(freq), dtype='complex')
    cdef double[:] ConvolvedMotion1 = np.empty(NFFT,dtype='float64')
    cdef double[:] ConvolvedMotion2 = np.empty(NFFT,dtype='float64')
    cdef double [:,:] RotConvolvedMotion = np.empty((N_rot_angle,LEN))
    cdef double[:,:] Sa = np.empty((len(T),len(RotD)),dtype='float64')
    cdef double[:] RotSaval = np.empty(N_rot_angle, dtype='float64')
    Fmotion1 = np.fft.rfft(np.asarray(motion1),NFFT)
    Fmotion2 = np.fft.rfft(np.asarray(motion2),NFFT)
    for k in range(len(T)):
        for u in range(len(freq)):
            FConvolvedMotion1[u] = Fmotion1[u]/(1.0+2.0j*D*freq[u]*T[k]-freq[u]*freq[u]*T[k]*T[k])
            FConvolvedMotion2[u] = Fmotion2[u]/(1.0+2.0j*D*freq[u]*T[k]-freq[u]*freq[u]*T[k]*T[k])
        ConvolvedMotion1 = np.fft.irfft(FConvolvedMotion1)
        ConvolvedMotion2 = np.fft.irfft(FConvolvedMotion2)
        RotConvolvedMotion = np.vstack((ConvolvedMotion1,ConvolvedMotion2)).T@rot
        RotSaval = np.max(RotConvolvedMotion, axis=0)
        for i in range(len(RotD)):
            Sa[k,i] = np.percentile(RotSaval, RotD[i])
    return(np.asarray(Sa))

def get_ngawest2_T():
    return(np.asarray([0.01, 0.02, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042, 0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.30, 0.32, 0.34, 0.35, 0.36, 0.38, 0.40, 0.42, 0.44, 0.45, 0.46, 0.48, 0.50, 0.55, 0.60, 0.65, 0.667, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0]))

def get_response_spectrum(*args, **kwargs):
    if('verbose' in kwargs):
        verbose = kwargs['verbose']
    else:
        verbose = 1
        
    if('T' in kwargs):
        T = np.asarray(kwargs['T'],dtype='float64')
    else:
        if(verbose != 0):
            print('using default T array from NGAwest2')
        T = get_ngawest2_T()
    if('motions' in kwargs):
        motions = np.asarray(kwargs['motions'],dtype='float64')
        if(len(motions.shape)==1):
            motions = motions.reshape((motions.size,1))
        if('RotD' in kwargs):
            if((motions.shape[0])!=2):
                print('You must specify two motions instead of ' + str(motions.shape[0]) + ' if you specify RotD')
                return
    else:
        print('You must specify motions')
        return
        
    if('D' in kwargs):
        D = kwargs['D']
    else:
        if(verbose != 0):
            print('setting D = 0.05')
        D=0.05
        
    if('zeropad' in kwargs):
        zeropad = kwargs['zeropad']
    else:
        if(verbose != 0):
            print('setting zeropad = 1')
        zeropad = 1

    if('dt' in kwargs):
        dt = kwargs['dt']
    else:
        print('You must specify dt')
        return
    
    if('RotD' in kwargs):
        RotD = np.asarray(kwargs['RotD'],dtype='float64')
        if(len(RotD)==1):
            RotD = RotD.reshape((RotD.size, 1))
        return(get_rotD_ARS(T, D, motions[0], motions[1], dt, zeropad, RotD, 100))
    else:
        return(get_ARS(T,D,motions,dt,zeropad))