import numpy as np
import ucla_geotech_tools.ipyconsol as pcl

# test 1, uniform soil constant load
def test_ipyconsol1():
    try:
        H = 4.35
        N = 100
        Cc = 5.729
        Cr = 0.532
        sigvref = 50.0
        esigvref = 8.763
        Gs = 1.686
        kref = 1.0e-7
        ekref = 8.344
        Ck = 200.0
        Ca = 0.0
        tref = 3070.0
        dsigv = 400.0
        ocrvoidratiotype = 0
        ocrvoidratio = 1.5
        ru = 0.0
        tmax = 35532000.0
        Ntime = 500
        qo = 10.0
        gammaw = 9.81
        pa = 101.325
        tol = 1.0e-8
        output = pcl.compute(
            Cc=Cc,
            Cr=Cr,
            sigvref=sigvref,
            esigvref=esigvref,
            Gs=Gs,
            kref=kref,
            ekref=ekref,
            Ck=Ck,
            Ca=Ca,
            tref=tref,
            qo=qo,
            dsigv=dsigv,
            ocrvoidratiotype=ocrvoidratiotype,
            ocrvoidratio=ocrvoidratio,
            drainagetype=0,
            ru=ru,
            pa=pa,
            tol=tol,
            N=N,
            H=H,
            Ntime=Ntime,
            tmax=tmax
            )
        print('test_ipyconsol.py test 1 finished')
    except:
        print('test_ipyconsol.py test 1 failed')

# test2 non-uniform soil, constant load
def test_ipyconsol2():
    try:
        H1 = 2.0
        H2 = 2.35
        N1 = 50
        N2 = 50
        depth1 = np.linspace(0,H1,num=N1,dtype=float,endpoint=False)
        depth2 = np.linspace(H1,H1+H2,num=N2+1,dtype=float)
        depth = np.hstack((depth1,depth2))
        N = depth.size-1
        Cc = np.full(N+1,5.729,dtype=float)
        Cr = np.full(N+1,0.532,dtype=float)
        sigvref = np.full(N+1,50,dtype=float)
        esigvref = np.full(N+1,8.763,dtype=float)
        Gs = np.full(N+1,1.868,dtype=float)
        kref = np.full(N+1,1.e-7,dtype=float)
        ekref = np.full(N+1,8.344,dtype=float)
        Ck = np.full(N+1,200,dtype=float)
        Ca = np.full(N+1,0.325,dtype=float)
        tref = np.full(N+1,3070,dtype=float)
        dsigv = np.full(N+1,27,dtype=float)
        ocrvoidratiotype = np.hstack((np.full(N1,0,dtype=int),np.full(N2+1,1,dtype=int)))
        ocrvoidratio = np.hstack((np.full(N1,1.5,dtype=float),np.full(N2+1,10.0,dtype=float)))
        ru = np.full(N+1,0,dtype=float)
        tmax = 35532000
        Ntime = 100
        qo=10
        gammaw = 9.81
        pa = 101.325
        tol = 1.0e-8

        OCR = []
        sigmap = []
        eo = []
        for i in range(depth.size):
            if(ocrvoidratiotype[i]==0):
                OCR.append(ocrvoidratio[i])
                sigmap.append(None)
                eo.append(None)
            if(ocrvoidratiotype[i]==1):
                eo.append(ocrvoidratio[i])
                OCR.append(None)
                sigmap.append(None)
            if(ocrvoidratiotype[i]==2):
                sigmap.append(ocrvoidratio[i])
                OCR.append(None)
                eo.append(None)

        depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype = pcl.get_inputs(Cc=Cc,Cr=Cr,sigvref=sigvref,esigvref=esigvref,Gs=Gs,kref=kref,ekref=ekref,Ck=Ck,Ca=Ca,tref=tref,qo=qo,dsigv=dsigv,ocrvoidratiotype=ocrvoidratiotype,ocrvoidratio=ocrvoidratio,drainagetype=0,ru=ru,pa=pa,tol=tol,depth=depth,Ntime=Ntime,tmax=tmax)
        output = pcl.compute(depth=depth,Cc=Cc,Cr=Cr,sigvref=sigvref,esigvref=esigvref,Gs=Gs,kref=kref,ekref=ekref,Ck=Ck,Ca=Ca,tref=tref,qo=qo,dsigv=dsigv,ocrvoidratiotype=ocrvoidratiotype,ocrvoidratio=ocrvoidratio,Ntime=Ntime,tmax=tmax,drainagetype=0,ru=ru,pa=pa,tol=tol)
        print('test_ipyconsol.py test 2 finished')
    except:
        print('test_ipyconsol.py test 2 failed')

# test 3 uniform soil, time-varying load
def test_ipyconsol3():
    try:
        H = 4.35
        N = 100
        Cc = 5.729
        Cr = 0.532
        sigvref = 50.0
        esigvref = 8.763
        Gs = 1.686
        kref = 1.0e-7
        ekref = 8.344
        Ck = 200.0
        Ca = 0.325
        tref = 3070.0
        dsigv = 27.0
        ocrvoidratiotype = 0
        ocrvoidratio = 1.5
        ru = 0.0
        tmax = 35532000.0
        Ntime = 100
        time = np.logspace(np.log10(tmax)-5,np.log10(tmax),num=Ntime)
        loadfactor = []
        for t in time:
            if(t<100000):
                loadfactor.append(t/100000)
            else:
                loadfactor.append(1.0)
        loadfactor = np.array(loadfactor, dtype=float)
        qo = 10.0
        gammaw = 9.81
        pa = 101.325
        tol = 1.0e-8

        OCR = []
        sigmap = []
        eo = []
        for i in range(N+1):
            if(ocrvoidratiotype==0):
                OCR.append(ocrvoidratio)
                sigmap.append(None)
                eo.append(None)
            if(ocrvoidratiotype==1):
                eo.append(ocrvoidratio)
                OCR.append(None)
                sigmap.append(None)
            if(ocrvoidratiotype==2):
                sigmap.append(ocrvoidratio)
                OCR.append(None)
                eo.append(None)
        output = pcl.compute(Cc=Cc,Cr=Cr,sigvref=sigvref,esigvref=esigvref,Gs=Gs,kref=kref,ekref=ekref,Ck=Ck,Ca=Ca,tref=tref,qo=qo,dsigv=dsigv,ocrvoidratiotype=ocrvoidratiotype,ocrvoidratio=ocrvoidratio,drainagetype=0,ru=ru,pa=pa,tol=tol,N=N,H=H,time=time,loadfactor=loadfactor)
        print('test_ipyconsol.py test 3 finished')
    except:
        print('test_ipyconsol.py test 3 failed')

# test 4 non-uniform soil, time-varying load
def test_ipyconsol4():
    try:
        H1 = 2.0
        H2 = 2.35
        N1 = 50
        N2 = 50
        Cc1 = 5.729
        depth1 = np.linspace(0,H1,num=N1,dtype=float,endpoint=False)
        depth2 = np.linspace(H1,H1+H2,num=N2+1,dtype=float)
        depth = np.hstack((depth1,depth2))
        N = depth.size-1
        Cc = np.full(N+1,5.729,dtype=float)
        Cr = np.full(N+1,0.532,dtype=float)
        sigvref = np.full(N+1,50,dtype=float)
        esigvref = np.full(N+1,8.763,dtype=float)
        Gs = np.full(N+1,1.868,dtype=float)
        kref = np.full(N+1,1.e-7,dtype=float)
        ekref = np.full(N+1,8.344,dtype=float)
        Ck = np.full(N+1,200,dtype=float)
        Ca = np.full(N+1,0.325,dtype=float)
        tref = np.full(N+1,3070,dtype=float)
        dsigv = np.full(N+1,27,dtype=float)
        ocrvoidratiotype = np.hstack((np.full(N1,0,dtype=int),np.full(N2+1,1,dtype=int)))
        ocrvoidratio = np.hstack((np.full(N1,1.5,dtype=float),np.full(N2+1,10.0,dtype=float)))
        ru = np.full(N+1,0,dtype=float)
        tmax = 35532000
        Ntime = 100
        time = np.logspace(np.log10(tmax)-5,np.log10(tmax),num=Ntime)
        loadfactor = []
        for t in time:
            if(t<100000):
                loadfactor.append(t/100000)
            else:
                loadfactor.append(1.0)
        loadfactor = np.array(loadfactor, dtype=float)
        qo=10
        gammaw = 9.81
        pa = 101.325
        tol = 1.0e-8

        OCR = []
        sigmap = []
        eo = []
        for i in range(depth.size):
            if(ocrvoidratiotype[i]==0):
                OCR.append(ocrvoidratio[i])
                sigmap.append(None)
                eo.append(None)
            if(ocrvoidratiotype[i]==1):
                eo.append(ocrvoidratio[i])
                OCR.append(None)
                sigmap.append(None)
            if(ocrvoidratiotype[i]==2):
                sigmap.append(ocrvoidratio[i])
                OCR.append(None)
                eo.append(None)
        inputs = [Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,dsigv,OCR,sigmap,eo]
        input_labels = [r'$C_c$',r'$C_r$',r'$\sigma_{v,ref}$',r'$e_{\sigma v,ref}$',r'$G_s$',r'$k_{ref}$',r'$e_{k,ref}$',r'$C_k$',r'$C_\alpha$',r'$t_{ref}$',r'$\Delta\sigma_v$',r'$OCR$',r"$\sigma_p'$",r'$e_o$']
                
        output = pcl.compute(depth=depth,Cc=Cc,Cr=Cr,sigvref=sigvref,esigvref=esigvref,Gs=Gs,kref=kref,ekref=ekref,Ck=Ck,Ca=Ca,tref=tref,qo=qo,dsigv=dsigv,ocrvoidratiotype=ocrvoidratiotype,ocrvoidratio=ocrvoidratio,time=time,loadfactor=loadfactor,drainagetype=0,ru=ru,pa=pa,tol=tol)
        print('test_ipyconsol.py test 4 finished')
    except:
        print('test_ipyconsol.py test 4 failed')