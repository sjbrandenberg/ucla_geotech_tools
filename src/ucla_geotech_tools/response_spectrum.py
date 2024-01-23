import numpy as np

def get_ngawest2_T():
    """Return numpy array containing spectral periods utilized in NGAWest2 project"""
    return(np.asarray([0.01, 0.02, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042, 0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.30, 0.32, 0.34, 0.35, 0.36, 0.38, 0.40, 0.42, 0.44, 0.45, 0.46, 0.48, 0.50, 0.55, 0.60, 0.65, 0.667, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0]))

def next_fast_len(target):
    """Return a number of datapoints for efficient calculation of fast Fourier transform. Used in zero padding"""
    I2 = 0
    while(2**I2<target):
        I2+=1
    I3 = 0
    while(3**I3<target):
        I3+=1
    I5 = 0
    while(5**I5<target):
        I5+=1
    return int(np.min([2**I2, 3**I3, 5**I5]))

def get_ARS(T,D,motions,dt,zeropad):
    """Return pseudo-acceleration response spectrum for a Numpy array containing acceleration records
    
    Inputs:
        T = numpy array of spectral periods in seconds
        D = damping for single-degree-of-freedom oscillator
        motions = Numpy array containing input acceleration records
        dt = time step in seconds
        zeropad = 0 for no zero padding, 1 for zero padding. Zero padding makes fast Fourier transform faster
    """
    if(zeropad==1):
        NFFT = next_fast_len(len(motions[0]))
    else:
        NFFT = len(motions[0])
    freq = np.fft.rfftfreq(NFFT,dt)
    Fmotions = np.fft.rfft(motions)
    Sa = [None]*len(motions)
    for i, Fmotion in enumerate(Fmotions):
        Fconvolved_motion = Fmotion/(1.0+2.0j*D*freq*T[:,np.newaxis]-freq**2*T[:,np.newaxis]**2)
        convolved_motion = np.fft.irfft(Fconvolved_motion)
        Sa[i] = np.max(np.absolute(convolved_motion),axis=1)
    return np.asarray(Sa)
                    
def get_rotD_ARS(T, D, motion1, motion2, dt, zeropad, RotD, N_rot_angle):
    """Return RotD pseudo-acceleration response spectrum for two Numpy arrays containing orthogonal horizontal acceleration records
    
    Inputs:
        T = numpy array of spectral periods in seconds
        D = damping for single-degree-of-freedom oscillator
        motion1 = Numpy array containing input acceleration record for motion 1
        motion2 = Numpy array containing input acceleration record for motion 2
        dt = time step in seconds
        zeropad = 0 for no zero padding, 1 for zero padding. Zero padding makes fast Fourier transform faster
        RotD = Numpy array containing desired RotD values (e.g., RotD = np.asarray([10,50,90]))
        N_rot_angle = Number of angles to use in computing RotD values
    """
    if(len(motion2) != len(motion1)):
        return('Both motions must have same length')
    if(zeropad==1):
        NFFT = next_fast_len(len(motion1))
    else:
        NFFT = len(motion1)
    theta = np.linspace(0,np.pi,N_rot_angle)
    rot = np.asarray([np.cos(theta), -np.sin(theta)])
    freq = np.fft.rfftfreq(NFFT,dt)
    Fmotion1 = np.fft.rfft(np.asarray(motion1),NFFT)
    Fmotion2 = np.fft.rfft(np.asarray(motion2),NFFT)
    
    for k in range(len(T)):
        FConvolvedMotion1 = Fmotion1/(1.0+2.0j*D*freq*T[:,np.newaxis]-freq**2*T[:,np.newaxis]**2)
        FConvolvedMotion2 = Fmotion2/(1.0+2.0j*D*freq*T[:,np.newaxis]-freq**2*T[:,np.newaxis]**2)
        ConvolvedMotion1 = np.fft.irfft(FConvolvedMotion1)
        ConvolvedMotion2 = np.fft.irfft(FConvolvedMotion2)
        RotConvolvedMotion = np.vstack((ConvolvedMotion1,ConvolvedMotion2)).T@rot
        RotSaval = np.max(RotConvolvedMotion, axis=0)
        for i in range(len(RotD)):
            Sa[k,i] = np.percentile(RotSaval, RotD[i])
    return(np.asarray(Sa))


def get_response_spectrum(**kwargs):
    """Function that reads keyword arguments and calls get_ARS function
    
    Keyword arugments:
        motions = numpy array containing acceleration time series (required)
        dt = time step for acceleration records (required)
        T = numpy array of spectral periods (default = NGAWest2 spectral periods)
        D = damping for single-degree-of-freedom oscillator 
        zeropad = 0 for no zero padding, 1 for zero padding
    """
    if('motions' not in kwargs):
        print('you must specify motions')
        return
    if('dt' not in kwargs):
        print('you must specify dt')
        return
    verbose = kwargs.get('verbose',1)
    T = kwargs.get('T',get_ngawest2_T())
    motions = kwargs['motions']
    if(len(motions.shape)==1):
        motions = motions.reshape((1,len(motions)))
    D = kwargs.get('D',0.05)
    zeropad = kwargs.get('zeropad',1)
    dt = kwargs['dt']
    N_rot_angle = kwargs.get('N_rot_angle',100)
    if('RotD' in kwargs):
        RotD = np.asarray(kwargs['RotD'],dtype=float)
        if(len(RotD)==1):
            RotD = RotD.reshape((1,RotD.size))
        return(get_rotD_ARS(T, D, motions[0], motions[1], dt, zeropad, RotD, N_rot_angle))
    else:
        return get_ARS(T,D,motions,dt,zeropad)
