import ucla_geotech_tools.auto_fchp as af
import numpy as np

dt = 0.02
acc = np.load('./tests/2D.OBS03..BH1__20110107T211622Z__20110107T212922Z.npy')

def test_auto_fchp1():
    try:
        fchp = af.get_fchp(dt=dt,acc=acc)
        print('test_fchp.py command 1 passed. fchp = ', fchp)
    except:
        print('test_fchp.py command 1 failed')
def test_auto_fchp2():
    try:
        fchp = af.get_fchp(dt=dt,acc=acc,target=0.02,tol=0.001,poly_order=6,maxiter=30,fchp_min=0.001,fchp_max=0.5,filter_order=5.0,tukey_alpha=0.05,apply_disp_ratio=1,disp_ratio_time=2,disp_ratio_target=0.02)
        print('test_fchp.py command 2 passed. fchp = ', fchp)
    except:
        print('test_fchp.py command 2 failed')