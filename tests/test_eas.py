import numpy as np
import pandas as pd
import ucla_geotech_tools.eas as eas

def test_eas():
    try:
        acc1 = pd.read_csv('20210102144223_BK_BDM_HHE.txt',sep=' ')['acc(g)'].values
        acc2 = pd.read_csv('20210102144223_BK_BDM_HHN.txt',sep=' ')['acc(g)'].values
        dt = 0.01
        b = 188.5
        w = 1.0/(10.0**(3.0/b))
        fc, eas_vec = eas.get_smooth_eas(acc1, acc2, dt, b, w)
        print('test_eas passed')
    except:
        print('test_eas failed')
