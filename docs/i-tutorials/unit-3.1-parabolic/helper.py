import scipy.sparse as sp
import matplotlib.pylab as plt
from ngsolve import *

def ShowPattern(A,precision=-1,binarize=False):
    rows_all,cols_all,vals_all = A.COO()
    for i,(r,c,v) in enumerate(zip(rows_all,cols_all,vals_all)):
        vals_all[i] = abs(vals_all[i])
        if binarize and vals_all[i] > precision:
            vals_all[i] = 1.0
    minval = 0
    maxval = max(vals_all)
    A_all = sp.csr_matrix((vals_all,(rows_all,cols_all)),shape=(A.height,A.width))
    plt.figure(figsize=(6,6*A.width/A.width))
    plt.imshow(A_all.todense(),interpolation='none',cmap='jet',vmin=minval,vmax=maxval)
    plt.show()
