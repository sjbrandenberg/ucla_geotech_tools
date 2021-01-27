import numpy as np
from libc.math cimport pow, log10, log, fabs, exp

def get_inputs(*args, **kwargs):
    cdef double[:] depth, Cc, Cr, sigvref, esigvref, Gs, kref, ekref, Ck, Ca, tref, dsigv, time, loadfactor, ru, ocrvoidratio
    cdef int[:] ocrvoidratiotype
    cdef double qo, pa, tol, gammaw, H
    cdef int N, drainagetype
    required = ['Cc','Cr','sigvref','esigvref','Gs','kref','ekref','Ck','Ca','tref','qo','dsigv','ocrvoidratiotype','ocrvoidratio']
    required_type= [np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.float, np.int, np.float]
    for req in required:
        if req not in kwargs:
            print('You must specify ' + req)
            return
    if(isinstance(kwargs['Cc'], (int,float))):
        for req in required:
            if(not isinstance(kwargs[req], (int,float))):
                print('If Cc is a number, ' + req + ' must also be a number')
                return
        required1 = ['N','H']
        for req1 in required1:
            if req1 not in kwargs:
                print('you must specify ' + req1 + ' when using input option 1, uniform soil')
                return
        N = int(kwargs['N'])
        H = float(kwargs['H'])
        Cc = np.full(N+1, kwargs['Cc'], dtype=np.float)
        Cr = np.full(N+1, kwargs['Cr'], dtype=np.float)
        sigvref = np.full(N+1, kwargs['sigvref'], dtype=np.float)
        esigvref = np.full(N+1, kwargs['esigvref'], dtype=np.float)
        Gs = np.full(N+1, kwargs['Gs'], dtype=np.float)
        kref = np.full(N+1, kwargs['kref'], dtype=np.float)
        ekref = np.full(N+1, kwargs['ekref'], dtype=np.float)
        Ck = np.full(N+1, kwargs['Ck'], dtype=np.float)
        Ca = np.full(N+1, kwargs['Ca'], dtype=np.float)
        tref = np.full(N+1, kwargs['tref'], dtype=np.float)
        qo = float(kwargs['qo'])
        dsigv = np.full(N+1, kwargs['dsigv'], dtype=np.float)
        ocrvoidratiotype = np.full(N+1, kwargs['ocrvoidratiotype'], dtype=np.int)
        ocrvoidratio = np.full(N+1, kwargs['ocrvoidratio'], dtype=np.float)
        depth = np.linspace(0,H,num=N+1,dtype=np.float,endpoint=True)
        if('ru' in kwargs):
            ru = np.full(N+1,kwargs['ru'],dtype=np.float)
        else:
            ru = np.full(N+1,0.0,dtype=np.float)
    elif(isinstance(kwargs['Cc'],(np.ndarray,(list,np.ndarray)))):
        for req in required:
            if req == 'qo':
                continue
            if len(kwargs[req]) != len(kwargs['Cc']):
                print('All vectors must have the same length using input option 2, layered soil. Cc has length = ' + str(len(kwargs['Cc'])), ' while ' + req + ' has length = ' + str(len(kwargs[req])))
                return
        Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio = [kwargs[req] for req in required]
        if('depth' in kwargs):
            depth = kwargs['depth']
        else:
            print('You must specify depth if you specify Cc as a list or array.')
            return
        if('ru' in kwargs):
            if(isinstance(kwargs['ru'],(list,np.ndarray))):
                ru = kwargs['ru']
            else:
                print('If Cc is a number, ru must also be a number')
        else:
            ru = np.full(kwargs['N']+1,0.0,dtype=np.float)
    else:
        print('Cc must be either a scalar, list, or np.ndarray')
        return
    if 'Ntime' in kwargs:
        if 'tmax' not in kwargs:
            print('You must specify tmax if you specify Ntime')
            return
        else:
            time = np.logspace(np.log10(kwargs['tmax'])-5,np.log10(kwargs['tmax']),num=kwargs['Ntime'],dtype=np.float)
            loadfactor = np.full(kwargs['Ntime'],1.0,dtype=np.float)
    elif 'time' in kwargs:
        if 'loadfactor' not in kwargs:
            print('You must specify loadfactor if you specify time')
            return
        if len(kwargs['loadfactor']) != len(kwargs['time']):
            print('time and loadfactor must have the same length. Currently time has length = ' + str(len(kwargs['time'])) + 'and loadfactor has length = ' + str(len(kwargs['loadfactor'])))
            return
        time = kwargs['time']
        loadfactor = kwargs['loadfactor']
    else:
        print('You must specify either Ntime and tmax, or time and loadfactor')
        return
    
    if 'gammaw' in kwargs:
        gammaw = kwargs['gammaw']
    else:
        gammaw = 9.81
    if 'tol' in kwargs:
        tol = kwargs['tol']
    else:
        tol = 1.0e-8
    if 'pa' in kwargs:
        pa = kwargs['pa']
    else:
        pa = 101.325
    if('drainagetype' in kwargs):
        drainagetype = kwargs['drainagetype']
    else:
        print('drainagetype not specified. Assuming double drainage.')
        drainagetype = 0
    return [depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype]
    
def Logspace(double tmin, double tmax, int Ntime, double tstart):
    cdef double dlogt = log10((tmax-tmin)/tstart)/(Ntime-1)
    cdef double logt = log10(tstart)
    cdef double[:] t = np.zeros(Ntime)
    for j in range(Ntime):
        t[j] = pow(10.0,logt) + tmin
        logt = dlogt + logt
    return t

def get_ktest(double etest,double ekref,double kref,double Ck):
    return kref*pow(10.0,(etest-ekref)/Ck)

def get_avtest(double elast, double sigvlast, double sigvtest, double eref, double sigvref, double Cc, double Cr):
    cdef double sigmap = pow(10.0,(eref-elast+Cc*log10(sigvref)-Cr*log10(sigvlast))/(Cc-Cr))
    cdef double De
    cdef double av
    if(sigmap<=sigvlast and sigvtest>sigvlast):
        De = Cc*log10(sigvtest/sigvlast)
    elif(sigmap>sigvtest or sigvtest<sigvlast):
        De = Cr*log10(sigvtest/sigvlast)
    else:
        De = Cr*log10(sigmap/sigvlast) + Cc*log10(sigvtest/sigmap)
    if(fabs((sigmap-sigvlast)/sigvlast)<1.e-5 and fabs((sigvlast-sigvtest)/sigvtest)<1.e-5):
        av = Cc/sigvlast/log(10.0)
    elif(sigmap>sigvtest and fabs((sigvlast-sigvtest)/sigvtest)<1.e-5):
        av = Cr/sigvlast/log(10.0)
    else:
        av = De/(sigvtest-sigvlast)
        print av
    return av

def get_etest(double etest, double elast, double eref, double alpha, double tref, double sigvtest, double sigvlast, double sigvref, double Cc, double avtest, double utest, double ulast, double dt):
    cdef double etol = 1.0e-12
    cdef int I
    cdef double res
    cdef double dres_detest
    if(alpha>0):
        I=0
        res = avtest*(utest-ulast) - 0.5*dt*alpha/tref*(exp((etest-eref)/alpha + Cc/alpha*log10(sigvtest/sigvref)) + exp((elast-eref)/alpha + Cc/alpha*log10(sigvlast/sigvref))) - etest + elast
        while(fabs(res)>etol):
            dres_detest = -0.5*dt/tref*exp((etest-eref)/alpha+Cc/alpha*log10(sigvtest/sigvref))-1.0
            etest = etest - res/dres_detest
            res = avtest*(utest-ulast) - 0.5*dt*alpha/tref*(exp((etest-eref)/alpha + Cc/alpha*log10(sigvtest/sigvref)) + exp((elast-eref)/alpha + Cc/alpha*log10(sigvlast/sigvref))) - etest + elast
            I+=1
            if(I==100):
                return -1
    else:
        etest = elast + avtest*(utest-ulast)
    return etest

def get_residual(double[:] ulast, double[:] utest, double[:] elast, double[:] etest, double[:] klast, double[:] ktest, double[:] avtest, double[:] zlast, double[:] ztest, double[:] sigvlast, double[:] sigvtest, double gammaw, double[:] Ca, double[:] Cc, double[:] sigvref, double[:] esigvref, double dt, double[:] tref, int drainagetype, int N):
    cdef double[:] residual = np.zeros(N+1)
    cdef double dzlast, dztest, res
    cdef int i
    for i in range(1,N,1):
        alpha = Ca[i]/log(10.0)
        dzlast = 0.5*(zlast[i+1]-zlast[i-1])
        dztest = 0.5*(ztest[i+1]-ztest[i-1])
        res = -avtest[i]/(1.0+elast[i])*(utest[i]-ulast[i])
        res += 0.5*dt/gammaw*(ktest[i]*(utest[i+1]-2.0*utest[i]+utest[i-1])/dztest/dztest+klast[i]*(ulast[i+1]-2.0*ulast[i]+ulast[i-1])/dzlast/dzlast)
        res += 0.125*dt/gammaw*((utest[i+1]-utest[i-1])*(ktest[i+1]-ktest[i-1])/dztest/dztest+(ulast[i+1]-ulast[i-1])*(klast[i+1]-klast[i-1])/dzlast/dzlast)
        if(alpha>0):
            res += 0.5*alpha*dt/tref[i]/(1+elast[i])*(exp((etest[i]-esigvref[i])/alpha+Cc[i]/alpha*log10(sigvtest[i]/sigvref[i])) + exp((elast[i]-esigvref[i])/alpha+Cc[i]/alpha*log10(sigvlast[i]/sigvref[i])))
        residual[i] = res

    if(drainagetype==0):
        residual[0] = 0.0
        residual[N] = 0.0

    if(drainagetype==1):
        residual[0] = 0.0
        dzlast = (zlast[N]-zlast[N-1])
        dztest = (ztest[N]-ztest[N-1])
        residual[N] = -avtest[N]/(1.0+elast[N])*(utest[N]-ulast[N])
        residual[N] += 0.5*dt/gammaw*(ktest[N]*(utest[N-1]-utest[N])/dztest/dztest+klast[N]*(ulast[N-1]-ulast[N])/dzlast/dzlast)
        residual[N] += 0.125*dt/gammaw*((utest[N]-utest[N-1])*(ktest[N]-ktest[N-1])/dztest/dztest+(ulast[N]-ulast[N-1])*(klast[N]-klast[N-1])/dzlast/dzlast)
        if(alpha>0):
            residual[N] += 0.5*alpha*dt/tref[i]/(1+elast[N])*(exp((etest[N]-esigvref[N])/alpha+Cc[i]/alpha*log10(sigvtest[N]/sigvref[i])) + exp((elast[N]-esigvref[N])/alpha+Cc[N]/alpha*log10(sigvlast[N]/sigvref[i])))

    if(drainagetype==2):
        residual[N] = 0.0
        dzlast = (zlast[1]-zlast[0])
        dztest = (ztest[1]-ztest[0])
        residual[0] = -avtest[0]/(1.0+elast[0])*(utest[0]-ulast[0])
        residual[0] += 0.5*dt/gammaw*(ktest[0]*(utest[1]-utest[0])/dztest/dztest+klast[0]*(ulast[1]-ulast[0])/dzlast/dzlast)
        residual[0] += 0.125*dt/gammaw*((utest[1]-utest[0])*(ktest[1]-ktest[0])/dztest/dztest+(ulast[1]-ulast[0])*(klast[1]-klast[0])/dzlast/dzlast)
        if(alpha>0):
            residual[0] += 0.5*alpha*dt/tref[i]/(1+elast[0])*(exp((etest[0]-esigvref[0])/alpha+Cc[i]/alpha*log10(sigvtest[0]/sigvref[0])) + exp((elast[0]-esigvref[0])/alpha+Cc[0]/alpha*log10(sigvlast[0]/sigvref[0])))
    return residual

def get_utest2(int N, double[:] klast, double[:] ktest, double[:] zlast, double[:] ztest, double[:] avtest, double[:] elast, double[:] etest, double[:] ulast, double[:] utest, double[:] sigvlast, double[:] sigvtest, double[:] Res, double dt, double[:] Cc, double[:] Cr, double[:] Ca, double gammaw, double[:] ekref, double[:] kref, double[:] Ck, double[:] tref, double[:] esigvref, double[:] sigvref, double[:] dsigv, int drainagetype, double pa):
    cdef int arrayLength, i
    cdef double dzlast, dztest, dresdutest, dresduim1test, dresduip1test, dresdavtest, dresdktest, dresdkim1test, dresdkip1test, dresdsigvtest, dresdetest, detestdutest, deim1testduim1test, deip1testduip1test, detestdsigvtest, dresddztest, dktestdetest, dkim1testdeim1test, dkip1testdeip1test, sigmap, davtestdsigvtest, ddztestdetest, dresdu, alpha
    if(drainagetype==0):
        arrayLength = N-1
    else:
        arrayLength = N
    cdef double[:] a = np.zeros(arrayLength)
    cdef double[:] b = np.zeros(arrayLength)
    cdef double[:] c = np.zeros(arrayLength)
    cdef double[:] x = np.zeros(arrayLength)
    cdef double[:] uout = np.zeros(len(utest))
    a[0] = 0.0
    c[arrayLength-1] = 0.0
    
    for i in range(1,N,1):
        alpha = Ca[i]/log(10.0)
        dzlast = 0.5*(zlast[i+1]-zlast[i-1])
        dztest = 0.5*(ztest[i+1]-ztest[i-1])
        dresdutest = -avtest[i]/(1.0+elast[i])-ktest[i]*dt/dztest/dztest/gammaw
        dresduim1test = 0.5*dt/gammaw/dztest/dztest*(ktest[i] + 0.25*(ktest[i-1] - ktest[i+1]))
        dresduip1test = 0.5*dt/gammaw/dztest/dztest*(ktest[i] + 0.25*(ktest[i+1] - ktest[i-1]))
        dresdavtest = (ulast[i]-utest[i])/(1.0+elast[i])
        dresdktest = 0.5*dt/gammaw/dztest/dztest*(utest[i-1]-2.0*utest[i]+utest[i+1])
        dresdkim1test = 0.125*dt/gammaw/dztest/dztest*(utest[i-1]-utest[i+1])
        dresdkip1test = -dresdkim1test
        if(alpha>0.0):
            dresdsigvtest = 0.5*dt*Cc[i]/log(10.0)/sigvtest[i]/tref[i]/(1.0+elast[i])*exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i]))
            dresdetest = 0.5*dt/tref[i]/(1.0+elast[i])*exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i]))
            detestdutest = avtest[i]/(1.0+0.5*dt/tref[i]*exp((etest[i]-esigvref[i])/alpha+Cc[i]/alpha*log10(sigvtest[i]/sigvref[i])))
            deim1testduim1test = avtest[i-1]/(1.0+0.5*dt/tref[i-1]*exp((etest[i-1]-esigvref[i-1])/alpha+Cc[i-1]/alpha*log10(sigvtest[i-1]/sigvref[i-1])))
            deip1testduip1test = avtest[i+1]/(1.0+0.5*dt/tref[i+1]*exp((etest[i+1]-esigvref[i+1])/alpha+Cc[i+1]/alpha*log10(sigvtest[i+1]/sigvref[i+1])))
            detestdsigvtest = 0.5*Cc[i]*dt/tref[i]/sigvtest[i]/log(10.0)*exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i]))/(1.0+0.5*dt/tref[i]*exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i])))
        else:
            dresdsigvtest = 0.0
            dresdetest = 0.0
            detestdutest = avtest[i]
            deim1testduim1test = avtest[i-1]
            deip1testduip1test = avtest[i+1]
            detestdsigvtest = 0.0
        dresddztest = -dt/gammaw/dztest/dztest/dztest*(ktest[i]*(utest[i-1]-2.0*utest[i]+utest[i+1]) + 0.25*(ktest[i+1]-ktest[i-1])*(utest[i+1]-utest[i-1]))
        dktestdetest = kref[i]*log(10.0)/Ck[i]*pow(10.0,(etest[i]-ekref[i])/Ck[i])
        dkim1testdeim1test = kref[i-1]*log(10.0)/Ck[i-1]*pow(10.0,(etest[i-1]-ekref[i-1])/Ck[i-1])
        dkip1testdeip1test = kref[i+1]*log(10.0)/Ck[i+1]*pow(10.0,(etest[i+1]-ekref[i+1])/Ck[i+1])
        sigmap = pow(10.0,(esigvref[i]-elast[i]+Cc[i]*log10(sigvref[i])-Cr[i]*log10(sigvlast[i]))/(Cc[i]-Cr[i]))
        
        if(fabs((sigvtest[i]-sigvlast[i])/sigvlast[i])>1.e-10*pa):
            if(sigmap<=sigvlast[i] and sigvtest[i] > sigvlast[i]):
                davtestdsigvtest = Cc[i]*(sigvtest[i] - sigvtest[i]*log(sigvtest[i]/sigvlast[i]) - sigvlast[i])/(sigvtest[i]*log(10.0)*(sigvtest[i]-sigvlast[i])*(sigvtest[i]-sigvlast[i]))
            elif(sigmap>sigvtest[i] or sigvtest[i] < sigvlast[i]):
                davtestdsigvtest = Cr[i]*(sigvtest[i] - sigvtest[i]*log(sigvtest[i]/sigvlast[i]) - sigvlast[i])/(sigvtest[i]*log(10.0)*(sigvtest[i]-sigvlast[i])*(sigvtest[i]-sigvlast[i]))
            else:
                davtestdsigvtest = (Cc[i]*(sigvtest[i] - sigvtest[i]*log(sigvtest[i]/sigmap) - sigvlast[i]) - Cr[i]*sigvtest[i]*log(sigmap/sigvlast[i]))/(sigvtest[i]*log(10.0)*(sigvtest[i]-sigvlast[i])*(sigvtest[i]-sigvlast[i]))
        else:
            if(fabs((sigvtest[i]-sigmap)/sigmap)<=1.e-10*pa):
                davtestdsigvtest = -Cc[i]*log(10.0)/sigvtest[i]/sigvtest[i]
            else:
                davtestdsigvtest = -Cr[i]*log(10.0)/sigvtest[i]/sigvtest[i]
        ddztestdetest = dzlast/(1.0+elast[i])
        dresdu = dresdutest - dresdsigvtest + dresdktest*dktestdetest*detestdutest - dresdavtest*davtestdsigvtest + dresddztest*ddztestdetest*detestdutest + dresdetest*(detestdutest - detestdsigvtest)
        if(drainagetype==0 or drainagetype==1):
            b[i-1] = dresdu
            if(i>1):
                a[i-1] = dresduim1test + dresdkim1test*dkim1testdeim1test*deim1testduim1test
            if(i<N-1):
                c[i-1] = dresduip1test + dresdkip1test*dkip1testdeip1test*deip1testduip1test
        if(drainagetype==2):
            b[i] = dresdu
            if(i>1):
                a[i] = dresduim1test + dresdkim1test*dkim1testdeim1test*deim1testduim1test
            if(i<N-1):
                c[i] = dresduip1test + dresdkip1test*dkip1testdeip1test*deip1testduip1test
    if(drainagetype==1):
        alpha = Ca[N]/log(10.0)
        dztest = 0.5*(ztest[N]-ztest[N-2])
        dzlast = 0.5*(zlast[N]-zlast[N-2])
        dresduip1test = 0.5*dt/gammaw/dztest/dztest*(ktest[N-1] + 0.25*(ktest[N]-ktest[N-2]))
        dresdkip1test = -0.5*dt/gammaw/dztest/dztest*(utest[N-2]-utest[N])
        dkip1testdeip1test = kref[N]*log(10.0)/Ck[N]*pow(10.0,(etest[N]-ekref[N])/Ck[N])
        if(alpha>0):
            deip1testduip1test = avtest[N]/(1.0+0.5*dt/tref[N]*exp((etest[N]-esigvref[N])/alpha+Cc[N]/alpha*log10(sigvtest[N]/sigvref[N])))
        else:
            deip1testduip1test = avtest[N]
        c[N-2] = dresduip1test + dresdkip1test*dkip1testdeip1test*deip1testduip1test
        dzlast = zlast[N] - zlast[N-1]
        dztest = ztest[N] - ztest[N-1]
        dresdutest = -avtest[N]/(1.0+elast[N]) - 0.125*dt/gammaw/dztest/dztest*(3.0*ktest[N]+ktest[N-1])
        dresduim1test = 0.125*dt/gammaw/dztest/dztest*(3.0*ktest[N]+ktest[N-1])
        dresdavtest = (ulast[N]-utest[N])/(1.0+elast[N])
        dresdktest = -0.375*dt/gammaw/dztest/dztest*(utest[N]-utest[N-1])
        dresdkim1test = -0.125*dt/gammaw/dztest/dztest*(utest[N]-utest[N-1])
        if(alpha>0.0):
            dresdsigvtest = 0.5*dt*Cc[N]/log(10.0)/sigvtest[N]/tref[N]/(1.0+elast[N])*exp((etest[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvtest[N]/sigvref[N]))
            dresdetest = 0.5*dt/tref[N]/(1.0+elast[N])*exp((etest[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvtest[N]/sigvref[N]))
            detestdutest = avtest[N]/(1.0+0.5*dt/tref[N]*exp((etest[N]-esigvref[N])/alpha+Cc[N]/alpha*log10(sigvtest[N]/sigvref[N])))
            deim1testduim1test = avtest[N-1]/(1.0+0.5*dt/tref[N]*exp((etest[N-1]-esigvref[N])/alpha+Cc[N]/alpha*log10(sigvtest[N-1]/sigvref[N])))
            deip1testduip1test = avtest[N]/(1.0+0.5*dt/tref[N]*exp((etest[N]-esigvref[N])/alpha+Cc[N]/alpha*log10(sigvtest[N]/sigvref[N])))
            detestdsigvtest = 0.5*Cc[N]*dt/tref[N]/sigvtest[N]/log(10.0)*exp((etest[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvtest[N]/sigvref[N]))/(1.0+0.5*dt/tref[N]*exp((etest[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvtest[N]/sigvref[N])))
        else:
            dresdsigvtest = 0.0
            dresdetest = 0.0
            detestdutest = avtest[N]
            deim1testduim1test = avtest[N-1]
            deip1testduip1test = avtest[N]
            detestdsigvtest = 0.0
        dresddztest = 0.25*dt/gammaw/dztest/dztest/dztest*(3.0*ktest[N]+ktest[N-1])*(utest[N]-utest[N-1])
        dktestdetest = kref[N]*log(10.0)/Ck[N]*pow(10.0,(etest[N]-ekref[N])/Ck[N])
        dkim1testdeim1test = kref[N-1]*log(10.0)/Ck[N-1]*pow(10.0,(etest[N-1]-ekref[N-1])/Ck[N-1])
        sigmap = pow(10.0,(esigvref[N]-elast[N]+Cc[N]*log10(sigvref[N])-Cr[N]*log10(sigvlast[N]))/(Cc[N]-Cr[N]))
        if(fabs((sigvtest[N]-sigvlast[N])/sigvlast[N])>1.e-10*pa):
            if(sigmap<=sigvlast[N] and sigvtest[N] > sigvlast[N]):
                davtestdsigvtest = Cc[N]*(sigvtest[N] - sigvtest[N]*log(sigvtest[N]/sigvlast[N]) - sigvlast[N])/(sigvtest[N]*log(10.0)*(sigvtest[N]-sigvlast[N])*(sigvtest[N]-sigvlast[N]))
            elif(sigmap>sigvtest[N] or sigvtest[N] < sigvlast[N]):
                davtestdsigvtest = Cr[N]*(sigvtest[N] - sigvtest[N]*log(sigvtest[N]/sigvlast[N]) - sigvlast[N])/(sigvtest[N]*log(10.0)*(sigvtest[N]-sigvlast[N])*(sigvtest[N]-sigvlast[N]))
            else:
                davtestdsigvtest = (Cc[N]*(sigvtest[N] - sigvtest[N]*log(sigvtest[N]/sigmap) - sigvlast[N]) - Cr[N]*sigvtest[N]*log(sigmap/sigvlast[N]))/(sigvtest[N]*log(10.0)*(sigvtest[N]-sigvlast[N])*(sigvtest[N]-sigvlast[N]))
        else:
            if(fabs((sigvtest[N]-sigmap)/sigmap)<=1.e-10*pa):
                davtestdsigvtest = -Cc[N]*log(10.0)/sigvtest[N]/sigvtest[N]
            else:
                davtestdsigvtest = -Cr[N]*log(10.0)/sigvtest[N]/sigvtest[N]
        ddztestdetest = dzlast/(1.0+elast[N])
        dresdu = dresdutest - dresdsigvtest + dresdktest*dktestdetest*detestdutest - dresdavtest*davtestdsigvtest + dresddztest*ddztestdetest*detestdutest + dresdetest*(detestdutest - detestdsigvtest)
        b[N-1] = dresdu
        a[N-1] = dresduim1test + dresdkim1test*dkim1testdeim1test*deim1testduim1test
    if(drainagetype==2):
        alpha = Ca[0]/log(10.0)
        dzlast = 0.5*(zlast[2] - zlast[0])
        dztest = 0.5*(ztest[2] - ztest[0])
        dresduim1test = 0.5*dt/gammaw/dztest/dztest*(ktest[1]+0.25*(ktest[0]-ktest[2]))
        dresdkim1test = 0.125*dt/gammaw/dztest/dztest*(utest[0]-utest[2])
        dkim1testdeim1test = kref[0]*log(10.0)/Ck[0]*pow(10.0,(etest[0]-ekref[0])/Ck[0])
        if(alpha>0):
            deim1testduim1test = avtest[0]/(1.0+0.5*dt/tref[0]*exp((etest[0]-esigvref[0])/alpha+Cc[0]/alpha*log10(sigvtest[0]/sigvref[0])))
        else:
            deim1testduim1test = avtest[0]
        a[1] = dresduim1test + dresdkim1test*dkim1testdeim1test*deim1testduim1test
        dzlast = zlast[1]-zlast[0]
        dztest = ztest[1]-ztest[0]
        dresdutest = -avtest[0]/(1.0+elast[0]) - 0.125*dt/gammaw/dztest/dztest*(3.0*ktest[0]+ktest[1])
        dresduip1test = 0.125*dt/gammaw/dztest/dztest*(3.0*ktest[0]+ktest[1])
        dresdktest = -0.375*dt/gammaw/dztest/dztest*(utest[0]-utest[1])
        dresdkip1test = -0.125*dt/gammaw/dztest/dztest*(ulast[0]-ulast[1])
        dresddztest = 0.25*dt/gammaw/dztest/dztest/dztest*(3.0*ktest[0]+ktest[1])*(utest[0]-utest[1])
        dresdavtest = (ulast[0]-utest[0])/(1.0+elast[0])
        if(alpha>0.0):
            dresdsigvtest = 0.5*dt*Cc[0]/log(10.0)/sigvtest[0]/tref[0]/(1.0+elast[0])*exp((etest[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvtest[0]/sigvref[0]))
            dresdetest = 0.5*dt/tref[0]/(1.0+elast[0])*exp((etest[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvtest[0]/sigvref[0]))
            detestdutest = avtest[0]/(1.0+0.5*dt/tref[0]*exp((etest[0]-esigvref[0])/alpha+Cc[0]/alpha*log10(sigvtest[0]/sigvref[0])))
            deim1testduim1test = avtest[0]/(1.0+0.5*dt/tref[0]*exp((etest[0]-esigvref[0])/alpha+Cc[0]/alpha*log10(sigvtest[0]/sigvref[0])))
            deip1testduip1test = avtest[1]/(1.0+0.5*dt/tref[1]*exp((etest[1]-esigvref[1])/alpha+Cc[1]/alpha*log10(sigvtest[1]/sigvref[1])))
            detestdsigvtest = 0.5*Cc[0]*dt/tref[0]/sigvtest[0]/log(10.0)*exp((etest[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvtest[0]/sigvref[0]))/(1.0+0.5*dt/tref[0]*exp((etest[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvtest[0]/sigvref[0])))
        else:
            dresdsigvtest = 0.0
            dresdetest = 0.0
            detestdutest = avtest[0]
            deim1testduim1test = avtest[0]
            deip1testduip1test = avtest[1]
            detestdsigvtest = 0.0
        dktestdetest = kref[0]*log(10.0)/Ck[0]*pow(10.0,(etest[0]-ekref[0])/Ck[0])
        dkip1testdeip1test = kref[1]*log(10.0)/Ck[1]*pow(10.0,(etest[1]-ekref[1])/Ck[1])
        sigmap = pow(10.0,(esigvref[0]-elast[0]+Cc[0]*log10(sigvref[0])-Cr[0]*log10(sigvlast[0]))/(Cc[0]-Cr[0]))
        if(fabs((sigvtest[0]-sigvlast[0])/sigvlast[0])>1.e-10*pa):
            if(sigmap<=sigvlast[0] and sigvtest[0] > sigvlast[0]):
                davtestdsigvtest = Cc[0]*(sigvtest[0] - sigvtest[0]*log(sigvtest[0]/sigvlast[0]) - sigvlast[0])/(sigvtest[0]*log(10.0)*(sigvtest[0]-sigvlast[0])*(sigvtest[0]-sigvlast[0]))
            elif(sigmap>sigvtest[0] or sigvtest[0] < sigvlast[0]):
                davtestdsigvtest = Cr[0]*(sigvtest[0] - sigvtest[0]*log(sigvtest[0]/sigvlast[0]) - sigvlast[0])/(sigvtest[0]*log(10.0)*(sigvtest[0]-sigvlast[0])*(sigvtest[0]-sigvlast[0]))
            else:
                davtestdsigvtest = (Cc[0]*(sigvtest[0] - sigvtest[0]*log(sigvtest[0]/sigmap) - sigvlast[0]) - Cr[0]*sigvtest[0]*log(sigmap/sigvlast[0]))/(sigvtest[0]*log(10.0)*(sigvtest[0]-sigvlast[0])*(sigvtest[0]-sigvlast[0]))
        else:
            if(fabs((sigvtest[0]-sigmap)/sigmap)<=1.e-10*pa):
                davtestdsigvtest = -Cc[0]*log(10.0)/sigvtest[0]/sigvtest[0]
            else:
                davtestdsigvtest = -Cr[0]*log(10.0)/sigvtest[0]/sigvtest[0]
        ddztestdetest = dzlast/(1.0+elast[0])
        dresdu = dresdutest - dresdsigvtest + dresdktest*dktestdetest*detestdutest - dresdavtest*davtestdsigvtest + dresddztest*ddztestdetest*detestdutest + dresdetest*(detestdutest - detestdsigvtest)
        b[0] = dresdu
        c[0] = dresduip1test + dresdkip1test*dkip1testdeip1test*deip1testduip1test
    if(drainagetype==0):
        for i in range(1,len(utest)-1,1):
            x[i-1] = Res[i]
    if(drainagetype==1):
        for i in range(1,len(utest),1):
            x[i-1] = Res[i]
    if(drainagetype==2):
        for i in range(0,len(utest)-1,1):
            x[i] = Res[i]
    c[0] = c[0]/b[0]
    x[0] = x[0]/b[0]
    for i in range(1,arrayLength,1):
        m = 1.0 / (b[i] - a[i]*c[i-1])
        c[i] = c[i]*m
        x[i] = (x[i] - a[i]*x[i-1])*m
    for i in range(arrayLength-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    if(drainagetype==0):
        for i in range(0,arrayLength,1):
            uout[i+1] = utest[i+1] - x[i]
        uout[0] = 0.0
        uout[N] = 0.0
    if(drainagetype==1):
        for i in range(0,arrayLength,1):
            uout[i+1] = utest[i+1] - x[i]
        uout[0] = 0.0
    if(drainagetype==2):
        for i in range(0,arrayLength,1):
            uout[i] = utest[i] - x[i]
        uout[N]=0.0
    return uout

def get_utest(double[:] klast, double[:] ktest, double[:] avtest, double[:] zlast, double[:] ztest, double[:] ulast, double[:] sigvlast, double[:] sigvtest, double[:] elast, double[:] etest, double[:] Ca, double[:] Cc, double[:] sigvref, double[:] esigvref, double dt, double[:] tref, int drainagetype, int N, double gammaw):
    cdef int arrayLength, i
    if(drainagetype==0):
        arrayLength = N-1
    elif(drainagetype==1 or drainagetype==2):
        arrayLength = N
    else:
        arrayLength = N+1
    cdef double[:] a = np.zeros(arrayLength)
    cdef double[:] b = np.zeros(arrayLength)
    cdef double[:] c = np.zeros(arrayLength)
    cdef double[:] x = np.zeros(arrayLength)
    cdef double[:] uout = np.zeros(N+1)
    cdef double alpha
    a[0] = 0.0
    c[arrayLength-1] = 0.0
    
    cdef double dzlast, dttest
    for i in range(1,N,1):
        alpha = Ca[i]/log(10.0)
        dzlast = 0.5*(zlast[i+1]-zlast[i-1])
        dztest = 0.5*(ztest[i+1]-ztest[i-1])
        if(i>1):
            a[i-1]=-0.5*dt/gammaw/dztest/dztest*(ktest[i]+0.25*(ktest[i-1]-ktest[i+1]))
        b[i-1]= avtest[i]/(1.0+elast[i]) + ktest[i]*dt/gammaw/dztest/dztest
        if(i<N-1):
            c[i-1]=-0.5*dt/gammaw/dztest/dztest*(ktest[i]-0.25*(ktest[i-1]-ktest[i+1]))
        if(alpha>0):
            x[i-1]=avtest[i]*ulast[i]/(1.0+elast[i]) + 0.5*klast[i]*dt/gammaw/dzlast/dzlast*(ulast[i-1]-2.0*ulast[i]+ulast[i+1])  + 0.5*alpha*dt/tref[i]/(1.0+elast[i])*(exp((elast[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvlast[i]/sigvref[i])) + exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i])))
            x[i-1]+= 0.125*dt/gammaw/dzlast/dzlast*(ulast[i+1]-ulast[i-1])*(klast[i+1]-klast[i-1])
        else:
            x[i-1]=avtest[i]*ulast[i]/(1.0+elast[i]) + 0.5*klast[i]*dt/gammaw/dzlast/dzlast*(ulast[i-1]-2.0*ulast[i]+ulast[i+1])
            x[i-1]+=0.125*dt/gammaw/dzlast/dzlast*(ulast[i+1]-ulast[i-1])*(klast[i+1]-klast[i-1])             
    if(drainagetype==1):
        alpha = Ca[N]/log(10.0)
        dztest = 0.5*(ztest[N] - ztest[N-2])
        dzlast = 0.5*(zlast[N] - zlast[N-2])
        c[N-2]= -0.5*dt/gammaw/dztest/dztest*(ktest[N-1]-0.25*(ktest[N-2]-ktest[N]))
        dztest = ztest[N] - ztest[N-1]
        dzlast = zlast[N] - zlast[N-1]
        b[N-1]=avtest[N]/(1.0+elast[N]) + 0.5*dt/gammaw*ktest[N]/dztest/dztest
        a[N-1]=-0.5*dt/gammaw/dztest/dztest*(ktest[N]+0.25*(ktest[N-1]-ktest[N]))
        if(alpha>0):
            x[N-1]=avtest[N]*ulast[N]/(1.0+elast[N]) + 0.5*klast[N]*dt/gammaw/dzlast/dzlast*(ulast[N-1]-2.0*ulast[N]+ulast[N])  + 0.5*alpha*dt/tref[N]/(1.0+elast[N])*(exp((elast[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvlast[N]/sigvref[N])) + exp((etest[N]-esigvref[N])/alpha + Cc[N]/alpha*log10(sigvtest[N]/sigvref[N])))
            x[N-1]+= 0.125*dt/gammaw/dzlast/dzlast*(ulast[N]-ulast[N-1])*(klast[N]-klast[N-1])
        else:
            x[N-1]=avtest[N]*ulast[N]/(1.0+elast[N]) + 0.5*klast[N]*dt/gammaw/dzlast/dzlast*(ulast[N-1]-2.0*ulast[N]+ulast[N])
            x[N-1]+= 0.125*dt/gammaw/dzlast/dzlast*(ulast[N]-ulast[N-1])*(klast[N]-klast[N-1])
            
    if(drainagetype==2):
        for i in range(1,N,1):
            alpha = Ca[i]/log(10.0)
            dzlast = 0.5*(zlast[i+1]-zlast[i-1])
            dztest = 0.5*(ztest[i+1]-ztest[i-1])
            if(i>1):
                a[i]=-0.5*dt/gammaw/dztest/dztest*(ktest[i]+0.25*(ktest[i-1]-ktest[i+1]))
            b[i]= avtest[i]/(1.0+elast[i]) + ktest[i]*dt/gammaw/dztest/dztest
            if(i<N-1):
                c[i]=-0.5*dt/gammaw/dztest/dztest*(ktest[i]+0.25*(ktest[i+1]-ktest[i-1]))
            if(alpha>0):
                x[i]=avtest[i]*ulast[i]/(1.0+elast[i]) + 0.5*klast[i]*dt/gammaw/dzlast/dzlast*(ulast[i-1]-2.0*ulast[i]+ulast[i+1])  + 0.5*alpha*dt/tref[i]/(1.0+elast[i])*(exp((elast[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvlast[i]/sigvref[i])) + exp((etest[i]-esigvref[i])/alpha + Cc[i]/alpha*log10(sigvtest[i]/sigvref[i])))
                x[i]+= 0.125*dt/gammaw/dzlast/dzlast*(ulast[i+1]-ulast[i-1])*(klast[i+1]-klast[i-1])
            else:
                x[i]=avtest[i]*ulast[i]/(1.0+elast[i]) + 0.5*klast[i]*dt/gammaw/dzlast/dzlast*(ulast[i-1]-2.0*ulast[i]+ulast[i+1])
                x[i]+= 0.125*dt/gammaw/dzlast/dzlast*(ulast[i+1]-ulast[i-1])*(klast[i+1]-klast[i-1])
        
        alpha = Ca[0]/log(10.0)
        dztest = 0.5*(ztest[2] - ztest[0])
        dzlast = 0.5*(zlast[2] - zlast[0])
        a[1]=-0.5*dt/gammaw/dztest/dztest*(ktest[1]+0.25*(ktest[2]-ktest[0]))
        dzlast = zlast[1] - zlast[0]
        dztest = ztest[1] - ztest[0]
        b[0]= avtest[0]/(1.0+elast[0]) + ktest[0]*dt/gammaw/dztest/dztest
        c[0]= -0.5*dt/gammaw/dztest/dztest*(ktest[0]+0.25*(ktest[1]-ktest[0]))
        if(alpha>0):
            x[0]=avtest[0]*ulast[0]/(1.0+elast[0]) + 0.5*dt/gammaw*(klast[1]-klast[0])*(ulast[1]-ulast[0])/dzlast/dzlast + 0.5*alpha*dt/tref[0]/(1.0+elast[0])*(exp((elast[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvlast[0]/sigvref[0])) + exp((etest[0]-esigvref[0])/alpha + Cc[0]/alpha*log10(sigvtest[0]/sigvref[0])))
            x[0]+=0.125*dt/gammaw/dzlast/dzlast*(ulast[1]-ulast[0])*(klast[1]-klast[0])
        else:
            x[0]=avtest[0]*ulast[0]/(1.0+elast[0]) + 0.5*dt/gammaw*(klast[1]-klast[0])*(ulast[1]-ulast[0])/dzlast/dzlast
            x[0]+=0.125*dt/gammaw/dzlast/dzlast*(ulast[1]-ulast[0])*(klast[1]-klast[0])
    
    c[0] = c[0]/b[0]
    x[0] = x[0]/b[0]
    for i in range(1,arrayLength,1):
        m = 1.0 / (b[i] - a[i]*c[i-1])
        c[i] = c[i]*m
        x[i] = (x[i] - a[i]*x[i-1])*m       

    for i in range(arrayLength-2,-1,-1):
        x[i] = x[i] - c[i]*x[i+1]
    
    if(drainagetype==0):
        uout[0] = 0.0
        uout[N] = 0.0
        for i in range(0,arrayLength,1):
            uout[i+1] = x[i]
    if(drainagetype==1):
        uout[0] = 0.0
        for i in range(0,arrayLength,1):
            uout[i+1] = x[i]
    if(drainagetype==2):
        uout[N] = 0.0
        for i in range(0,arrayLength,1):
            uout[i] = x[i]
    return uout

def Max(double[:] v):
    return np.max([np.max(v),-np.min(v)])

def get_initial(depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype):
    cdef double[:] sigvo = np.zeros(depth.size)
    cdef double[:] sigvo2 = np.zeros(depth.size)
    cdef double[:] sigvf = np.zeros(depth.size)
    cdef double[:] eo = np.zeros(depth.size)
    cdef double[:] sigvncl = np.zeros(depth.size)
    sigvo[0] = qo
    sigvo2[0] = qo*(1-ru[0])
    sigvf[0] = qo+dsigv[0]*loadfactor[0] 
    if(ocrvoidratiotype[0]==0):
        eo[0] = esigvref[0]-Cc[0]*log10(ocrvoidratio[0]*qo/sigvref[0]) + Cr[0]*log10(ocrvoidratio[0])
        sigvncl[0] = pow(10.0,(esigvref[0] - eo[0] + Cc[0]*log10(sigvref[0]) - Cr[0]*log10(sigvo[0]))/(Cc[0]-Cr[0]))
    elif(ocrvoidratiotype[0]==1):
        eo[0] = ocrvoidratio[0]
        sigvncl[0] = pow(10.0,(esigvref[0] - eo[0] + Cc[0]*log10(sigvref[0]) - Cr[0]*log10(sigvo[0]))/(Cc[0]-Cr[0]))
    elif(ocrvoidratiotype[0]==2):
        eo[0] = esigvref[0] - Cc[0]*log10(ocrvoidratio[0]/sigvref[0]) + Cr[0]*log10(ocrvoidratio[0]/sigvo[0])
        sigvncl[0] = pow(10.0,(esigvref[0] - eo[0] + Cc[0]*log10(sigvref[0]) - Cr[0]*log10(sigvo[0]))/(Cc[0]-Cr[0]))
    else:
        print('Invalid ocrvoidratiotype: ' + str(ocrvoidratiotype[0]))
        return
    
    for i in range(1,depth.size,1):
        dz = depth[i] - depth[i-1]
        sigvotest = sigvo[i-1] + (Gs[i]-1.0)/(1.0+eo[i-1])*gammaw*dz
        if(ocrvoidratiotype[i]==0):
            eotest = esigvref[i] - Cc[i]*log10(ocrvoidratio[i]*sigvotest/sigvref[i])+Cr[i]*log10(ocrvoidratio[i])
            res=sigvotest-sigvo[i-1]-(Gs[i-1]-1)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
            while(fabs(res)>1.0e-10*pa):
                sigvotest = sigvo[i-1] + (Gs[i]-1.0)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
                eotest = esigvref[i] - Cc[i]*log10(ocrvoidratio[i]*sigvotest/sigvref[i])+Cr[i]*log10(ocrvoidratio[i])
                res=sigvotest-sigvo[i-1]-(Gs[i]-1)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
            sigvo[i] = sigvotest
            eo[i] = eotest
            sigvo2[i] = sigvo[i]*(1-ru[i])
            sigvf[i] = sigvo[i]+dsigv[i]*loadfactor[0]
            sigvncl[i] = pow(10.0,(esigvref[i] - eo[i] + Cc[i]*log10(sigvref[i]) - Cr[i]*log10(sigvo[i]))/(Cc[i]-Cr[i]))
        elif(ocrvoidratiotype[i]==1):
            eo[i] = ocrvoidratio[i]
            sigvo[i] = sigvo[i-1] + (Gs[i]-1.0)*dz*gammaw/(1.0+0.5*(eo[i]+eo[i-1]))
            sigvo2[i] = sigvo[i]*(1-ru[i])
            sigvf[i] = sigvo[i]+dsigv[i]*loadfactor[0]
            sigvncl[i] = pow(10.0,(esigvref[i] - eo[i] + Cc[i]*log10(sigvref[i]) - Cr[i]*log10(sigvo[i]))/(Cc[i]-Cr[i]))
        elif(ocrvoidratiotype[i]==2):
            sigvotest = sigvo[i-1] + (Gs[i]-1.0)/(1.0+eo[i-1])*gammaw*dz
            eotest = esigvref[i] - Cc[i]*log10(ocrvoidratio[i]/sigvref[i])+Cr[i]*log10(ocrvoidratio[i]/sigvotest)
            res=sigvotest-sigvo[i-1]-(Gs[i]-1)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
            while(fabs(res)>1.0e-10*pa):
                sigvotest = sigvo[i-1] + (Gs[i]-1.0)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
                eotest = esigvref[i] - Cc[i]*log10(ocrvoidratio[i]/sigvref[i])+Cr[i]*log10(ocrvoidratio[i]/sigvotest)
                res=sigvotest-sigvo[i-1]-(Gs[i]-1)*gammaw*dz/(1.0+0.5*(eotest+eo[i-1]))
            sigvo[i] = sigvotest
            eo[i] = eotest
            sigvo2[i] = sigvo[i]*(1-ru[i])
            sigvf[i] = sigvo[i]+dsigv[i]*loadfactor[0]
            sigvncl[i] = pow(10.0,(esigvref[i] - eo[i] + Cc[i]*log10(sigvref[i]) - Cr[i]*log10(sigvo[i]))/(Cc[i]-Cr[i]))
        else:
            print('Invalid ocrvoidratiotype: ' + str(ocrvoidratiotype[0]))
            return
    return {'eo':eo, 'sigvo':sigvo, 'sigvo2':sigvo2, 'sigvf':sigvf}
    #return([eo,sigvo,sigvo2,sigvf,sigvncl])

 
def compute(*args, **kwargs):
    cdef double[:] depth, Cc, Cr, sigvref, esigvref, Gs, kref, ekref, Ck, Ca, tref, dsigv, time, loadfactor, ru, ocrvoidratio 
    cdef int[:] ocrvoidratiotype
    cdef double qo, pa, tol, gammaw, H
    cdef int N, Ntime, drainagetype, i, j
    depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype = get_inputs(*args, **kwargs)
    
    N = depth.size - 1
    Ntime = time.size
    H = depth[N]
             
    cdef double[:] dsigvsub = np.zeros(N+1, dtype=float)
    cdef double[:] dsigvsublast = np.zeros(N+1, dtype=float)
    cdef double[:,:] u = np.zeros((N+1,Ntime))
    cdef double[:,:] sigv = np.zeros((N+1,Ntime))
    cdef double[:,:] z = np.zeros((N+1,Ntime))
    cdef double[:,:] e = np.zeros((N+1,Ntime))
    cdef double[:] sigvo = np.zeros(N+1)
    cdef double[:] sigvo2 = np.zeros(N+1)
    cdef double[:] sigvf = np.zeros(N+1)
    cdef double[:] eo = np.zeros(N+1)
    cdef double[:] sigvncl = np.zeros(N+1)
    cdef double[:] utest = np.zeros(N+1)
    cdef double[:] ulast = np.zeros(N+1)
    cdef double[:] sigvtest = np.zeros(N+1)
    cdef double[:] sigvlast = np.zeros(N+1)
    cdef double[:] avtest = np.zeros(N+1)
    cdef double[:] etest = np.zeros(N+1)
    cdef double[:] elast = np.zeros(N+1)
    cdef double[:] ktest = np.zeros(N+1)
    cdef double[:] klast = np.zeros(N+1)
    cdef double[:] epsvtest = np.zeros(N+1)
    cdef double[:] epsvlast = np.zeros(N+1)
    cdef double[:] ztest = np.zeros(N+1)
    cdef double[:] zlast = np.zeros(N+1)
    cdef double[:] Res = np.zeros(N+1)
    cdef double sigvotest, eotest, res, alpha, Hdr, cv, Tmax, Umax
    
    initial = get_initial(depth,Cc,Cr,sigvref,esigvref,Gs,kref,ekref,Ck,Ca,tref,qo,dsigv,ocrvoidratiotype,ocrvoidratio,ru,time,loadfactor,gammaw,tol,pa,drainagetype)
    eo = initial['eo']
    sigvo = initial['sigvo']
    sigvo2 = initial['sigvo2']
    sigvf = initial['sigvf']
            
    for i in range(0,N+1,1):
        dsigvsub[i] = dsigv[i]*loadfactor[0]
        dsigvsublast[i] = dsigv[i]*loadfactor[0]
        alpha = Ca[i]/log(10)
        ulast[i] = sigvf[i]-sigvo2[i]
        sigvlast[i] = sigvf[i] - ulast[i]
        avtest[i] = get_avtest(elast[i],sigvo2[i],sigvlast[i],esigvref[i],sigvref[i],Cc[i],Cr[i])
        elast[i] = eo[i] - avtest[i]*(sigvlast[i]-sigvo2[i])
        klast[i] = get_ktest(elast[i],ekref[i],kref[i],Ck[i])
    zlast[N] = H
    for i in range(0,N,1):
        zlast[N-1-i] = zlast[N-i]-H/N
    for j in range(0,Ntime,1):
        for i in range(0,N+1,1):
            dsigvsub[i] = dsigv[i]*loadfactor[j]
            ulast[i] = ulast[i] + dsigvsub[i] - dsigvsublast[i]
            dsigvsublast[i] = dsigvsub[i]
            ktest[i] = klast[i]
            ztest[i] = zlast[i]
            sigvtest[i] = sigvlast[i]
            etest[i] = elast[i]
            utest[i] = ulast[i]
            avtest[i] = get_avtest(elast[i],sigvlast[i],sigvtest[i],esigvref[i],sigvref[i],Cc[i],Cr[i])
        if(j==0):
            dt = time[j]
        else:
            dt = time[j] - time[j-1]
        utest = get_utest(klast,ktest,avtest,zlast,ztest,ulast,sigvlast,sigvtest,elast,etest,Ca,Cc,sigvref,esigvref,dt,tref,drainagetype,N,gammaw)
        for i in range(0,N+1,1):
            alpha = Ca[i]/log(10.0)
            sigvtest[i] = sigvlast[i] + ulast[i] - utest[i]
            avtest[i] = get_avtest(elast[i], sigvlast[i], sigvtest[i], esigvref[i], sigvref[i], Cc[i], Cr[i])
            if(alpha>0):
                etest[i] = get_etest(etest[i], elast[i], esigvref[i], alpha, tref[i], sigvtest[i], sigvlast[i], sigvref[i], Cc[i], avtest[i], utest[i], ulast[i], dt)
                epsvtest[i] = (elast[i] - etest[i])/(1.0+elast[i])
            else:
                etest[i] = elast[i] - avtest[i]*(ulast[i]-utest[i])
                epsvtest[i] = (elast[i] - etest[i])/(1.0+elast[i])
            ktest[i] = get_ktest(etest[i], ekref[i], kref[i], Ck[i])
        ztest[N] = H
        for i in range(0,N,1):
            ztest[N-1-i] = ztest[N-i] - (zlast[N-i] - zlast[N-i-1] - (epsvtest[N-i-1])*(zlast[N-i]-zlast[N-i-1]))
        Res = get_residual(ulast, utest, elast, etest, klast, ktest, avtest, zlast, ztest, sigvlast, sigvtest, gammaw, Ca, Cc, sigvref, esigvref, dt, tref, drainagetype,N)
        if(Max(Res)<=tol):
            for i in range(0,len(utest),1):
                u[i][j]=utest[i]
                ulast[i] = utest[i]
                sigvlast[i] = sigvtest[i]
                elast[i] = etest[i]
                klast[i] = ktest[i]
                zlast[i] = ztest[i]
                z[i][j]=ztest[i]
                sigv[i][j]=sigvtest[i]
                e[i][j] = etest[i]
            continue
        I=0
        while(Max(Res)>tol):
            utest=get_utest2(N,klast, ktest, zlast, ztest, avtest, elast, etest, ulast, utest, sigvlast, sigvtest, Res, dt, Cc, Cr, Ca, gammaw, ekref, kref, Ck, tref, esigvref, sigvref, dsigvsub, drainagetype, pa)
            for i in range(0,N+1,1):
                sigvtest[i] = sigvlast[i] + ulast[i] - utest[i]
                avtest[i] = get_avtest(elast[i], sigvlast[i], sigvtest[i], esigvref[i], sigvref[i], Cc[i], Cr[i])
                etest[i] = get_etest(etest[i], elast[i], esigvref[i], alpha, tref[i], sigvtest[i], sigvlast[i], sigvref[i], Cc[i], avtest[i], utest[i], ulast[i], dt)
                epsvtest[i] = (elast[i] - etest[i])/(1.0+elast[i])
                ktest[i] = get_ktest(etest[i], ekref[i], kref[i], Ck[i])
            
            ztest[N] = H
            for i in range(0,N,1):
                ztest[N-1-i] = ztest[N-i] - (zlast[N-i] - zlast[N-i-1] - (epsvtest[N-i-1])*(zlast[N-i]-zlast[N-i-1]))
            Res = get_residual(ulast, utest, elast, etest, klast, ktest, avtest, zlast, ztest, sigvlast, sigvtest, gammaw, Ca, Cc, sigvref, esigvref, dt, tref, drainagetype,N)
            if(Max(Res)<=tol):
                for i in range(0,len(utest),1):                    
                    u[i][j]=utest[i]
                    ulast[i] = utest[i]
                    sigvlast[i] = sigvtest[i]
                    elast[i] = etest[i]
                    klast[i] = ktest[i]
                    zlast[i] = ztest[i]
                    z[i][j]=ztest[i]
                    sigv[i][j]=sigvtest[i]  
                    e[i][j]=etest[i]
                continue         
            I=I+1
            print(Max(Res))
            if(I==100):
                print("ERROR! Increment " + str(j) + " did not converge in " + str(I) + " iterations.")
                return {'z':z, 'u':u, 'sigv':sigv, 'e':e}
    return {'z':z, 'u':u, 'sigv':sigv, 'e':e}
    #return([z,u,sigv,e,time,depth,eo,sigvo,sigvf])
