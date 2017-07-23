#Last-modified: 28 Apr 2014 02:31:22

from javelin.spear_covfunc import spear_covfunc as SCF
# from javelin.spear_covfunc import spear_covfunc as SCF_single
#from spear_covfunc import spear_covfunc as SCF
import py_spear_covfunc as py_SCF

import numpy as np

from javelin.threadpool import get_threadpool_size, map_noreturn
from javelin.gp import isotropic_cov_funs
from javelin.gp.GPutils import regularize_array

import unittest

""" The SPEAR covariance function, wrapper for the Fortran version.
"""


def spear_threading(x,y,idx,idy,sigma,tau,lags,wids,scales,scale_hidden=None,symm=None,set_pmap=False,blocksize=10000,baldwin=False) :
    """
    threaded version, divide matrix into subblocks with *blocksize*
    elements each. Do not use it when multiprocessing is on (e.g., in emcee MCMC
    sampling).

    set_pmap needs to be turned on for the Pmap_Model.
    """
    if (sigma<0. or tau<0.) :
        raise ValueError, 'The amp and scale parameters must be positive.'
    if (symm is None) :
        symm = (x is y) and (idx is idy)
    scale_hidden = np.array([scales[2 * k] for k in range(1, (len(scales) - 1) / 2 + 1)])

    x = regularize_array(x)
    y = regularize_array(y)
    nx = x.shape[0]
    ny = y.shape[0]
    if np.isscalar(idx) :
        idx = np.ones(nx, dtype="int", order="F")*idx
    if np.isscalar(idy) :
        idy = np.ones(nx, dtype="int", order="F")*idy
    # Figure out how to divide job up between threads (along y)
    n_threads = min(get_threadpool_size(), nx*ny/blocksize)
    if n_threads > 1 :
        if not symm:
            # divide ny evenly if x is not y
            bounds = np.linspace(0,ny,n_threads+1)
        else :
            # divide ny*ny evenly in quadrature if x is y
            bounds = np.array(np.sqrt(np.linspace(0,ny*ny,n_threads+1)),dtype=int)
    # Allocate the matrix
    C = np.asmatrix(np.empty((nx,ny),dtype=float,order='F'))
    if set_pmap :
        def targ(C,x,y,idx,idy,cmin,cmax,symm) :
            SCF.covmatpmap_bit(C,x,y,idx,idy,sigma,tau,lags,wids,scales,scale_hidden,1,cmin,cmax,symm)
    else :
        def targ(C,x,y,idx,idy,cmin,cmax,symm) :
            SCF.covmat_bit(C,x,y,idx,idy,sigma,tau,lags,wids,scales,cmin,cmax,symm)
    if n_threads <= 1 :
        targ(C,x,y,idx,idy,0,-1,symm)
    else :
        thread_args = [(C,x,y,idx,idy,bounds[i],bounds[i+1],symm) for i in xrange(n_threads)]
        map_noreturn(targ, thread_args)
    if symm:
        isotropic_cov_funs.symmetrize(C)
    # print "C: ", C
    return(C)


def spear(x,y,idx,idy,sigma,tau,lags,wids,scales,scale_hidden=None,symm=None,set_pmap=False,baldwin=False) :
    """ Clean version without multithreading. Used when multiprocessing is on
    (e.g., in emcee MCMC sampling).

    set_pmap needs to be turned on for the Pmap_Model.
    """
    if (sigma<0. or tau<0.) :
        raise ValueError, 'The amp and scale parameters must be positive.'
    if (symm is None) :
        symm = (x is y) and (idx is idy)
    
    # print "scale_hidden: ", scale_hidden
    

    # lags = np.array([0.0, 15.0, 0.0])
    # wids = np.array([0.0, 2.0, 0.0])
    # scales = np.array([1.0, 0.05, 0.9])
    if baldwin == False:
        scale_hidden = np.array([scales[2 * k] for k in range(1, (len(scales) - 1) / 2 + 1)])
    else:
        # print "scale_hidden: ", scale_hidden
        # print "scales: ", scales
        if scale_hidden is None:
            raise ValueError, 'unsupported scale_hidden type.'

    # sigma = 0.15
    # tau = 20.0
 
    # print "x: ", x
    # print "y: ", y
    # print "idx: ", idx
    # print "idy: ", idy


    x = regularize_array(x)
    y = regularize_array(y)
    nx = x.shape[0]
    ny = y.shape[0]
    if np.isscalar(idx) :
        idx = np.ones(nx, dtype="int", order="F")*idx
    if np.isscalar(idy) :
        idy = np.ones(nx, dtype="int", order="F")*idy
    # Allocate the matrix
    C = np.asmatrix(np.empty((nx,ny),dtype=float,order='F'))
    C_single = np.asmatrix(np.empty((nx,ny),dtype=float,order='F'))
    if set_pmap :
        if baldwin == True:
            SCF.covmatpmap_bit_baldwin(C,x,y,idx,idy,sigma,tau,lags,wids,scales,scale_hidden,1,0,-1,symm)
        else:
            SCF.covmatpmap_bit(C,x,y,idx,idy,sigma,tau,lags,wids,scales,scale_hidden,1,0,-1,symm)
    else :
        SCF.covmat_bit(C,x,y,idx,idy,sigma,tau,lags,wids,scales,0,-1,symm)
    if symm:
        isotropic_cov_funs.symmetrize(C)
        # isotropic_cov_funs.symmetrize(C_single)
    # print "C == C_single? ", (C == C_single).all()
    # yyy = input()
    return(C)



def spear_threading2(x,y,idx,idy,A,gamma,lags,wids,scales,scale_hidden=None,symm=None,set_pmap=False,blocksize=10000,baldwin=False) :
    """
    threaded version, divide matrix into subblocks with *blocksize*
    elements each. Do not use it when multiprocessing is on (e.g., in emcee MCMC
    sampling).

    set_pmap needs to be turned on for the Pmap_Model.
    """
    if (A<0. or gamma<0.) :
        raise ValueError, 'The amp and scale parameters must be positive.'
    if (symm is None) :
        symm = (x is y) and (idx is idy)
    x = regularize_array(x)
    y = regularize_array(y)
    nx = x.shape[0]
    ny = y.shape[0]
    if np.isscalar(idx) :
        idx = np.ones(nx, dtype="int", order="F")*idx
    if np.isscalar(idy) :
        idy = np.ones(nx, dtype="int", order="F")*idy

    scale_hidden = np.array([scales[2 * k] for k in range(1, (len(scales) - 1) / 2 + 1)])

    # Figure out how to divide job up between threads (along y)
    n_threads = min(get_threadpool_size(), nx*ny/blocksize)
    if n_threads > 1 :
        if not symm:
            # divide ny evenly if x is not y
            bounds = np.linspace(0,ny,n_threads+1)
        else :
            # divide ny*ny evenly in quadrature if x is y
            bounds = np.array(np.sqrt(np.linspace(0,ny*ny,n_threads+1)),dtype=int)
    # Allocate the matrix
    C = np.asmatrix(np.empty((nx,ny),dtype=float,order='F'))
    if set_pmap :
        def targ(C,x,y,idx,idy,cmin,cmax,symm) :
            SCF.covmatpmap_bit(C,x,y,idx,idy,A,gamma,lags,wids,scales,scale_hidden,2,cmin,cmax,symm)
    else :
        def targ(C,x,y,idx,idy,cmin,cmax,symm) :
            SCF.covmat_bit(C,x,y,idx,idy,A,gamma,lags,wids,scales,cmin,cmax,symm)
    if n_threads <= 1 :
        targ(C,x,y,idx,idy,0,-1,symm)
    else :
        thread_args = [(C,x,y,idx,idy,bounds[i],bounds[i+1],symm) for i in xrange(n_threads)]
        map_noreturn(targ, thread_args)
    if symm:
        isotropic_cov_funs.symmetrize(C)
    return(C)


def spear2(x,y,idx,idy,A,gamma,lags,wids,scales,scale_hidden=None,symm=None,set_pmap=False,baldwin=False) :
    """ Clean version without multithreading. Used when multiprocessing is on
    (e.g., in emcee MCMC sampling).

    set_pmap needs to be turned on for the Pmap_Model.
    """
    if (A<0. or gamma<0.) :
        raise ValueError, 'The amp and scale parameters must be positive.'
    if (symm is None) :
        symm = (x is y) and (idx is idy)
    x = regularize_array(x)
    y = regularize_array(y)
    nx = x.shape[0]
    ny = y.shape[0]
    if np.isscalar(idx) :
        idx = np.ones(nx, dtype="int", order="F")*idx
    if np.isscalar(idy) :
        idy = np.ones(nx, dtype="int", order="F")*idy
    
    scale_hidden = np.array([scales[2 * k] for k in range(1, (len(scales) - 1) / 2 + 1)])

    # Allocate the matrix
    C = np.asmatrix(np.empty((nx,ny),dtype=float,order='F'))
    if set_pmap :
        if baldwin == True:
            SCF.covmatpmap_bit_baldwin(C,x,y,idx,idy,A,gamma,lags,wids,scales,scale_hidden,2,0,-1,symm)
        else:
            SCF.covmatpmap_bit(C,x,y,idx,idy,A,gamma,lags,wids,scales,scale_hidden,2,0,-1,symm)
    else :
        SCF.covmat_bit2(C,x,y,idx,idy,A,gamma,lags,wids,scales,0,-1,symm)
    if symm:
        isotropic_cov_funs.symmetrize(C)
    # print "C: ", C
    return(C)





#I didn't modify this class. @ZHW
class PmapCovTest(unittest.TestCase):
    def testPmapCov(self):
        # fixed
        jdarr = np.array([0, 1, 0, 1]) # not sorted here.
        idarr = np.array([1, 1, 2, 2])
        # parameters, can be changed
        tau   = 1.0
        sigma = 1.0
        lags  = np.array([0.00, 0.25, 0.00])
        wids  = np.array([0.00, 0.00, 0.00])
        scales= np.array([1.00, 3.00, 2.00])
        C_true= np.empty((4,4), order="F")
        # diagonal
        C_true[0, 0] = 1.0
        C_true[1, 1] = 1.0
        C_true[2, 2] = scales[2]**2 + scales[2]*scales[1]*np.exp(-lags[1]/tau) + scales[1]*scales[2]*np.exp(-lags[1]/tau) + scales[1]**2
        C_true[3, 3] = scales[2]**2 + scales[2]*scales[1]*np.exp(-lags[1]/tau) + scales[1]*scales[2]*np.exp(-lags[1]/tau) + scales[1]**2
        # off
        C_true[0, 1] = C_true[1, 0] = np.exp(-1/tau)
        C_true[0, 2] = C_true[2, 0] = scales[2] + scales[1]*np.exp(-lags[1]/tau)
        C_true[0, 3] = C_true[3, 0] = scales[2]*np.exp(-1/tau) + scales[1]*np.exp(-(1.0-lags[1])/tau)
        C_true[1, 2] = C_true[2, 1] = scales[2]*np.exp(-1/tau) + scales[1]*np.exp(-(1.0+lags[1])/tau)
        C_true[1, 3] = C_true[3, 1] = scales[2] + scales[1]*np.exp(-lags[1]/tau)
        C_true[2, 3] = C_true[3, 2] = scales[2]*scales[2]*np.exp(-1/tau) + scales[2]*scales[1]*np.exp(-(1.0-lags[1])/tau) + scales[2]*scales[1]*np.exp(-(1.0+lags[1])/tau) + scales[1]*scales[1]*np.exp(-1/tau)
        C_true = C_true * sigma * sigma
        print "Truth :"
        print C_true
        # calculate from spear
        C_thread = spear_threading(jdarr, jdarr, idarr, idarr, sigma, tau, lags, wids, scales, symm=None, set_pmap=True)
        print "Calculated (threading) :"
        print C_thread
        C_bare = spear(jdarr, jdarr, idarr, idarr, sigma, tau, lags, wids, scales, symm=None, set_pmap=True)
        print "Calculated (no threading) :"
        print C_bare
        # compare
        self.assertTrue(np.allclose(C_true, C_thread, rtol=1e-05, atol=1e-08))
        self.assertTrue(np.allclose(C_true, C_bare,   rtol=1e-05, atol=1e-08))

if __name__ == "__main__":
    npt = 5
    ids = [1, 2]
    jd1 = np.linspace(0, 10, npt)
    jd2 = jd1[:]
    id1 = np.random.choice(ids, size=npt)
    id2 = id1[:]
    line_ind = np.where(id1 == 1)[0]
    cont_ind = np.where(id2 == 2)[0]
    line_mag = 22.3
    line_sigma = 2.0 
    cont_mag = 22.4
    cont_sigma = 1.0
    mag = np.zeros(npt)
    slagarr = np.array([0.0, 10.0, 0.0])
    swidarr = np.array([0.0, 2.0, 0.0])
    scalearr = np.array([1.0, 0.05, 0.89])
    scale_hidden = np.array([0.89])
    sigma = 0.15
    tau = 20.0
    C = np.asmatrix(np.empty((npt,npt),dtype=float,order='F'))
    for i in range(npt):
        if i in line_ind:
            mag[i] = np.random.normal(loc=line_mag, scale=line_sigma)
        elif i in cont_ind:
            mag[i] = np.random.normal(loc=cont_mag, scale=cont_sigma)


    print "indices are: ", id1
    # python covmatpmap result
    py_SCF.covmatpmap_bit(C, jd1, jd2, id1, id2, sigma, tau, slagarr, swidarr, scalearr,len(jd1),len(jd2),\
                          0, -1, True)
    print "the result of the python code is:\n ", C
    C = np.asmatrix(np.empty((npt,npt),dtype=float,order='F'))
    print "scalearr: ", scalearr
    f_SCF.covmatpmap_bit(C, jd1, jd2, id1, id2, sigma, tau, slagarr, swidarr, scalearr,scale_hidden,1,\
                          0, -1, True)
    isotropic_cov_funs.symmetrize(C)
    print "the result of the fortran code is:\n ", C

