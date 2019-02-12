# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:37:14 2017

@author: olem
"""

import scipy as sci
import numpy as np
from numba import jit
from numba import prange
from pytictoc import TicToc
import tables, warnings
from scipy import sparse
import h5py

def oem(y, Kin, xa, Seinvin, Sainvin, maxiter):
    K = Kin.todense()
    Seinv = Seinvin.todense()
    Sainv = Sainvin.todense()
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    xhat = sci.linalg.solve(M,b) + xa
                                   
    return xhat

def oem_basic(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    xhat = sci.sparse.linalg.spsolve(M,b) + xa
                                   
    return xhat

def oem_basic_sparse(y, K, xa, Seinv, Sainv, maxiter):
    
    xhat = sci.sparse.linalg.spsolve(((K.T.dot(Seinv)).dot(K) + Sainv).tocsc(),(K.T.dot(Seinv)).dot(y-K.dot(xa)))
    xhat = xhat + xa[:,0]                           
    return xhat

def oem_basic_sparse_2(y, K, xa, Seinv, Sainv, maxiter):
    
    S = (K.T.dot(Seinv)).dot(K) + Sainv
    KSe = (K.T).dot(Seinv)
    KSey = KSe.dot(y-K.dot(xa))
    xhat = sci.sparse.linalg.spsolve(S,KSey)                              
    xhat = xhat + xa[:,0]
    
    return xhat


def oem_cg(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x_0 = np.zeros(xa.shape)
    xhat = sci.sparse.linalg.cg(M,b,x0=x_0,maxiter=maxiter)[0] + xa[:,0]
                                   
    return xhat

def oem_cgs(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x_0 = np.zeros(xa.shape)
    xhat = sci.sparse.linalg.cgs(M,b,x0=x_0,maxiter=maxiter)[0] + xa
                                   
    return xhat

def oem_bicgstab(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x_0 = np.zeros(xa.shape)
    xhat = sci.sparse.linalg.bicgstab(M,b,x0=x_0,maxiter=maxiter)[0] + xa
                                   
    return xhat

def oem_gmres(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x_0 = np.zeros(xa.shape)
    xhat = sci.sparse.linalg.gmres(M,b,x0=x_0,maxiter=maxiter)[0] + xa
                                   
    return xhat

def oem_minres(y, K, xa, Seinv, Sainv, maxiter):
    
    M = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x_0 = np.zeros(xa.shape)
    xhat = sci.sparse.linalg.minres(M,b,x0=x_0,maxiter=maxiter)[0] + xa
                                   
    return xhat


def oem_cg2(y, K, xa, Seinv, Sainv,maxiter):  
    
    A = (K.T.dot(Seinv)).dot(K) + Sainv   
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x = np.zeros(xa.shape)
#    Ax = (K.T.dot(Seinv).dot(K) + Sainv).dot(x)
    Ax = A.dot(x)
    r = b-Ax
    
    p = r     
    rsold = r.T.dot(r)
    
    for i in range(1,maxiter):
        #Ap = K.T.dot((Seinv.dot(K.dot(p)))) + Sainv.dot(p)
        #t.tic()
        Ap = A.dot(p)
        #times[0] = t.tocvalue()
        #t.tic()
        #Ap = csrMult_parallel(A,p)
        #times[1] = t.tocvalue()
        #t.tic()
        alpha = rsold/(p.T.dot(Ap))
        #times[2] = t.tocvalue()
        #t.tic()
        x = x + alpha*p
        #times[3] = t.tocvalue()
        #t.tic()
        r = r - alpha*Ap
        #times[4] = t.tocvalue()
        #t.tic()
        rsnew = r.T.dot(r)
        #times[5] = t.tocvalue()
        #t.tic()
        p = (rsnew/rsold)*p + r
        #times[6] = t.tocvalue()
        #t.tic()
        rsold=rsnew
        #times[7] = t.tocvalue()
    
    xhat = x + xa
    
    return xhat

def oem_cg2_p(y, K, xa, Seinv, Sainv,maxiter):

    #K = Kin
    #xa = xain
    #Sainv = Sainvin
    
    
    
    A = (K.T.dot(Seinv)).dot(K) + Sainv
    b = ((K.T).dot(Seinv)).dot(y-K.dot(xa))
    x = np.zeros(xa.shape)
#    Ax = (K.T.dot(Seinv).dot(K) + Sainv).dot(x)
    Ax = A.dot(x)
    r = b-Ax
    
    p = r     
    rsold = r.T.dot(r)
    
    for i in range(1,maxiter):
        #Ap = K.T.dot((Seinv.dot(K.dot(p)))) + Sainv.dot(p)
        #t.tic()
        #Ap = A.dot(p)
        #times[0] = t.tocvalue()
        #t.tic()
        Ap = calc_Ax_p(Sainv,Seinv,K,p)
        #times[1] = t.tocvalue()
        #t.tic()
        alpha = rsold/(p.T.dot(Ap))
        #times[2] = t.tocvalue()
        #t.tic()
        x = x + alpha*p
        #times[3] = t.tocvalue()
        #t.tic()
        r = r - alpha*Ap
        #times[4] = t.tocvalue()
        #t.tic()
        rsnew = r.T.dot(r)
        #times[5] = t.tocvalue()
        #t.tic()
        p = (rsnew/rsold)*p + r
        #times[6] = t.tocvalue()
        #t.tic()
        rsold=rsnew
        #times[7] = t.tocvalue()
    
    xhat = x + xa
    
    return xhat

def oem_cg3(y, K, xa, Seinv, Sainv,maxiter):
    t = TicToc()
    Kt = K.T
    b = (Kt @ Seinv) @ (y-K @ xa)
    x = np.zeros(xa.shape)
    #Ax = (K.T.dot(Seinv).dot(K) + Sainv).dot(x)
    #Ax = A.dot(x)
    Ax = calc_Ax_2(Sainv,Seinv,K,Kt,x)
    r = b-Ax
    
    p = r     
    rsold = r.T.dot(r)
    for i in range(1,maxiter):
        #Ap = K.T.dot((Seinv.dot(K.dot(p)))) + Sainv.dot(p)
        t.tic()
        Ap =  calc_Ax_2(Sainv,Seinv,K,Kt,p)
        #Ap = csrMult_parallel(A,p)
        #times[1] = t.tocvalue()
        #t.tic()
        alpha = rsold/(p.T.dot(Ap))
        #print('alpha = ' + str(t.tocvalue()))
        #times[2] = t.tocvalue()
        #t.tic()
        x = x + alpha*p
        #print(t.tocvalue())
        #times[3] = t.tocvalue()
        #t.tic()
        r = r - alpha*Ap
        #print(t.tocvalue())
        #times[4] = t.tocvalue()
        #t.tic()
        rsnew = r.T.dot(r)
        #print(t.tocvalue())
        #times[5] = t.tocvalue()
        #t.tic()
        p = (rsnew/rsold)*p + r
        #print(t.tocvalue())
        #times[6] = t.tocvalue()
        #t.tic()
        rsold=rsnew
        #print(t.tocvalue())
        #times[7] = t.tocvalue()
    
    xhat = x + xa
    
    return xhat

def cg2(b,A,x0,maxiter):
    x = x0
    b = b
    Ax = A.dot(x)
    r = b-Ax
    
    p = r     
    rsold = r.T.dot(r)
    
    for i in range(1,maxiter):
        Ap = A.dot(p)
        alpha = rsold/(p.T.dot(Ap))
        x = x + alpha*p
        r = r - alpha*Ap
        rsnew = r.T.dot(r)
        p = (rsnew/rsold)*p + r
        rsold=rsnew
        
    xhat = x
    
    return xhat

def calc_Ax(Sainv,Seinv,K,x):    
    Mx = Sainv.dot(x) + (K.transpose()).dot(Seinv.dot(K.dot(x)));
    return Mx

def calc_Ax_2(Sainv,Seinv,K,Kt,x):
	t = TicToc()
    
	t.tic()
	Sainvx = Sainv @ x
	print('Sainvx ' + str(t.tocvalue()))
	
	t.tic()
	Kx = (K @ x)
	print('Kx ' + str(t.tocvalue()))
	
	t.tic()
	SeinvKx = Seinv @ Kx
	print('SeinvKx ' + str(t.tocvalue()))
	
	t.tic()
	KtSeinvKx = Kt @ SeinvKx
	print('KtSeinvKx ' +str(t.tocvalue()))

	t.tic()
	Mx = Sainvx + KtSeinvKx
	print('Sainv + KtSeinvKx ' + str(t.tocvalue()))
	return Mx

def calc_Ax_p(Sainv,Seinv,K,x):    
    Mx = csrMult_parallel(Sainv,x) + csrMult_parallel(K.transpose(),(csrMult_parallel(Seinv,(csrMult_parallel(K,x)))))
    Mx=np.expand_dims(Mx, axis=1)
    return Mx

def csrMult_parallel(A,x): 

    Adata= A.data
    Aindices= A.indices
    Aindptr=A.indptr
    Ashape=A.shape

    Ax = csrMult_parallel_execute(A.data,A.indices,A.indptr,A.shape,x)

    return Ax

@jit(parallel=True,nopython=True)
def csrMult_parallel_execute(Adata,Aindices,Aindptr,Ashape,x): 

    numRowsA = Ashape[0]    
    Ax = np.zeros(numRowsA)

    for i in prange(numRowsA):
        Ax_i = np.array([0.0])
        for dataIdx in range(Aindptr[i],Aindptr[i+1]):

            j = Aindices[dataIdx]
            Ax_i += Adata[dataIdx]*x[j,:]

        Ax[i] = Ax_i[0]            

    return Ax

#@jit(nopython=True, parallel=True)
def csrMult_parallel_time(A,x): 
    t = TicToc()    
    t.tic()
    Adata= A.data
    Aindices= A.indices
    Aindptr=A.indptr
    Ashape=A.shape
    times = t.tocvalue()
    numRowsA = Ashape[0]    
    Ax = np.zeros(numRowsA)

    for i in prange(numRowsA):
        Ax_i = 0.0        
        for dataIdx in range(Aindptr[i],Aindptr[i+1]):

            j = Aindices[dataIdx]
            Ax_i += Adata[dataIdx]*x[j]

        Ax[i] = Ax_i

    return Ax,times

def load_sparse_matrix(matname,size,fname='matlab.mat') :
    f = h5py.File(fname, 'r')
    data = f['/' + matname + '/data']
    jc = f['/' + matname + '/jc']
    ir = f['/' + matname + '/ir']
    M = sparse.csc_matrix((data.value, ir.value.astype(int), jc.value.astype(int)),shape=size,dtype='float64')
    M = M.tocsr()
    f.close()
    return M

def load_matrix(fname,matname) :
    f = h5py.File(fname, 'r')
    M = f[matname].value
    f.close()
    return M

#def tikonov(x,y,z,a0,ax,ay,az)   :
#    #Performs tikonov regularization on a matrix Sa
#    
#    nx = x.size
#    ny = y.size
#    nz = z.size
#    dx = x[1]-x[0]
#    dy = y[1]-y[0]
#    dz = z[1]-z[0]
#    
#    [L0, Lx, Ly, Lz] = calc_L(nx,ny,nz,dx,dy,dz)
#    
#    Sainv = a0[0].^2*(L0.T).dot(L0)
#    
#    return Sainv
#
#def calc_L(nx,ny,nz,dx,dy,dz)
#    
#    #Do not fill last row in Sa (due to numerical derivative)
#    
#    iP_all = zeros((nx-1)*(ny-1)*(nz-1),1);
#    iP_x_all = zeros((nx-1)*(ny-1)*(nz-1)*2,1);
#    iP_y_all = zeros((nx-1)*(ny-1)*(nz-1)*2,1);
#    iP_z_all = zeros((nx-1)*(ny-1)*(nz-1)*2,1);
#    
#    iP_derivative = zeros((nx-1)*(ny-1)*(nz-1)*2,1);
#    
#    %define grid positions
#    n=1;
#    m=1;
#    for ix=1:nx-1
#        %disp([num2str(ix) ' of ' num2str(nx)])
#        for iy=1:ny-1
#            for iz = 1:nz-1
#                iP = ix + (iy-1)*nx + (iz-1)*nx*ny; %find position of elenent
#                iP_x = (ix+1) + (iy-1)*nx + (iz-1)*nx*ny; %find position of first derivative
#                iP_y = ix + ((iy+1)-1)*nx + (iz-1)*nx*ny; %find position of first derivative
#                iP_z = ix + (iy-1)*nx + ((iz+1)-1)*nx*ny; %find position of first derivative
#                
#                iP_all(n) = iP;
#                iP_x_all(m) = iP;
#                iP_y_all(m) = iP;
#                iP_z_all(m) = iP;
#                iP_derivative(m) = -1;
#                m=m+1;
#                iP_x_all(m) = iP_x;
#                iP_y_all(m) = iP_y;
#                iP_z_all(m) = iP_z;
#                iP_derivative(m) = 1;
#                n=n+1;
#                m=m+1;
#            end
#        end
#    end
#    
#    index = [iP_all'; iP_all'];
#    
#    L0 = sparse(iP_all,iP_all,ones(size(iP_all)),nx*ny*nz,nx*ny*nz);
#    Lx = sparse(index(:),iP_x_all,iP_derivative,nx*ny*nz,nx*ny*nz);
#    Ly = sparse(index(:),iP_y_all,iP_derivative,nx*ny*nz,nx*ny*nz);
#    Lz = sparse(index(:),iP_z_all,iP_derivative,nx*ny*nz,nx*ny*nz);
#    
#    Lx = 1/dx*Lx;
#    Ly = 1/dy*Ly;
#    Lz = 1/dz*Lz;
#    
#    return L0, Lx, Ly, Lz
