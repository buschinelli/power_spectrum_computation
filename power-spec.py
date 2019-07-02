import numpy as np
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt
import random as rand

def power_spec(box,dx,bins):
    N = int(np.sqrt(np.size(box)))
    ibox = (dx**2)*fftshift(fft2(box))
    ibox_abs_sq = np.abs(ibox)**2
    K = k_matrix(box,dx)
    unique_K = np.unique(K,return_counts=True)
    max_k = np.amax(unique_K[0])
    ps = np.zeros((bins))
    step = max_k/bins
    p = 0
    while p*step < max_k:
        n_ele = 0
        for i in range(-N//2,N//2):
            for j in range(-N//2,N//2):
                z = np.sqrt(i**2+j**2)/(N*dx)
                if z < (p+1)*step and z >= p*step :
                    ps[p] = ps[p] + ibox_abs_sq[i,j]     
                    n_ele = n_ele + 1         
        ps[p] = ps[p]/n_ele
        p = p + 1
    
    ps = ps/((N*dx)**2)
    return ps


def k_matrix(box,dx):
    size = int(np.sqrt(np.size(box)))
    kx = fftshift(fftfreq(size,d=dx))
    ky = fftshift(fftfreq(size,d=dx))
    k = np.zeros((size,size),dtype='float')
    for i in range(0,size):
        for j in range(0,size):
            k[i,j] = np.sqrt(kx[i]**2+ky[j]**2)
    return k


def find_box(ps,dx,N):
    bins = len(ps)
    a = np.random.normal(size=(N,N))
    b = np.random.normal(size=(N,N))
    box = a + 1j*b
    shift_box = fftshift(box)
    K = k_matrix(box,dx)
    unique_K = np.unique(K,return_counts=True) 
    max_k = np.amax(unique_K[0])
    step = max_k/bins
    p = 0
    for i in range(-N//2,N//2):
        for j in range(-N//2,N//2):
            n_ele = 0
            while p*step < max_k:
                z = np.sqrt(i**2+j**2)/(N*dx)
                if z < (p+1)*step and z >= p*step:
                    el = np.random.normal(0,np.sqrt(ps[p]/2),size=(1,1))
                    box[i,j] = el[0]
                p = p + 1

    
    for i in range(0,N):
        for j in range(0,N):
            if box[i,j] != box[-i,-j]:
                box[i,j] = np.conjugate(box[-i,-j])
            else:
                box[i,j] = np.real(box[i,j])
                
    
    return ifft2(box)/(dx**2)
    
dx = 1
N = 4
bins = 4

box = np.random.normal(size=(N,N))
ps = power_spec(box,dx,bins)
boxf = find_box(ps,dx,N)
psf = power_spec(boxf,dx,bins)
plt.plot(make_k(box,dx,bins),ps)
plt.plot(make_k(box,dx,bins),psf)
