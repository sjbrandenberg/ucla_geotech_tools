import numpy as np
import scipy as sp
import ucla_geotech_tools.response_spectrum as rs

def test_response_spectrum1():
    try:
        N = 2**14
        dt = 0.005
        Nmotion = 20
        time = np.linspace(0,dt*N,N)
        f = np.linspace(0.5,8.0,Nmotion)
        motions = np.sin(2*np.pi*f[:,np.newaxis]*time)
        window = sp.signal.windows.tukey(N,0.20)
        mean_motions = np.average(motions,weights=window,axis=1)
        motions = window*(motions - mean_motions[:,np.newaxis])
        Fmotions = np.fft.rfft(motions)
        freq = np.fft.rfftfreq(N,dt)
        T = rs.get_ngawest2_T()
        D = 0.05
        zeropad = 0
        Sa = rs.get_response_spectrum(motions=motions,D=D,T=T,dt=dt,zeropad=0)
        print('test_response_spectrum1 passed')
    except:
        print('test_response_spectrum1 failed')

def test_response_spectrum2():
    try:
        N = 2**14
        dt = 0.005
        Nmotion = 2
        time = np.linspace(0,dt*N,N)
        f = np.linspace(0.5,8.0,Nmotion)
        motions = np.sin(2*np.pi*f[:,np.newaxis]*time)
        window = sp.signal.windows.tukey(N,0.20)
        mean_motions = np.average(motions,weights=window,axis=1)
        motions = window*(motions - mean_motions[:,np.newaxis])
        Fmotions = np.fft.rfft(motions)
        freq = np.fft.rfftfreq(N,dt)
        T = rs.get_ngawest2_T()
        D = 0.05
        zeropad = 0
        Sa = rs.get_response_spectrum(motions=motions,D=D,T=T,dt=dt,zeropad=1,RotD=50,N_rot_angle=100)
        Sa = rs.get_response_spectrum(motions=motions,D=D,T=T,dt=dt,zeropad=0,RotD=[0,50,100],N_rot_angle=100)
        print('test_response_spectrum1 passed')
    except:
        print('test_response_spectrum1 failed')
