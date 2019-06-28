import numpy as np
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt
import random as rand




## Take a NxN box in real space, compute its FFT and take its absolute
## value.

def abs_field(box):
    fft_box = fftshift(fft2((box)))
    psi_sqr = (np.abs(fft_box)**2)

    return psi_sqr


## Creates a grid containing all points in fourier space, and the values
## of k^2 = kx^2 + ky^2 corresponding to that point

def k_matrix(box):
    size = int(m.sqrt(np.size(box)))
    kx = fftshift(fftfreq(size))
    ky = fftshift(fftfreq(size))
    k = np.zeros((size,size),dtype='float')
    for i in range(0,size):
        for j in range(0,size):
            k[i,j] = np.sqrt(kx[i]**2+ky[j]**2)
    return k

## Outputs all possible values for k in a array, from smallar to larger k's

def unique_k(k):
    return np.unique(k)

## Computes the 1D power spectrum array, containing P(k), for each value of k

def power_spec(box):
    k = k_matrix(box)
    size = int(m.sqrt(np.size(box)))
    unique_elements = np.unique(k)
    counts_elements = np.unique(k, return_counts=True)
    fft_psi_abs = abs_field(box)
    kx = fftshift(fftfreq(size))
    ky = fftshift(fftfreq(size))
    ps = np.zeros(len(unique_elements),dtype='float')
    
    for i in range(0,len(unique_elements)):
        for j in range(0,len(kx)):
            for n in range(0,len(ky)):
                if unique_elements[i] == m.sqrt(kx[j]**2 + ky[n]**2):
                    ps[i] = ps[i] + fft_psi_abs[j,n]
        ps[i] = ps[i]/counts_elements[1][i]    
    return ps

## Plots the power spectrum as a function of k.

def plot_ps(box):
    ps = power_spec(box)
    k = unique_k(k_matrix(box))
    plt.plot(k,ps)
    plt.xlabel('Frequency')
    plt.ylabel('Power spectrum')
    plt.show()
    

############### REVERSE PROCESS #################
#################################################

        



def count_dim(ps_size):
    for i in range(1,1001):
        if ps_size == count_size(2*i):
            n = int(2*i)
            return n

    return 0


## This is the main code: it uses as input a given power spectrum and compute some box that would have as power spec
## the array ps inputed that way

def find_box(ps):
    size = np.size(ps)
    N = count_dim(size)
    a = np.random.normal(size=(N,N))
    b = np.random.normal(size=(N,N))
    box = a + 1j*b
    shift_box = fftshift(box)
    K = fftshift(k_matrix(box))
    k_value = np.unique(K,return_counts=True)    

    
    for i in range(0,N):
        for j in range(0,N):
            k = np.sqrt(i**2+j**2)/N
            k_compare = k_value[0] - k
            loc_k = np.where(k_compare == 0)[0]
            if len(loc_k) > 0 :
                box[i,j] = np.sqrt(ps[loc_k[0]]/2)*box[i,j]
    
    for i in range(0,N):
        for j in range(0,N):
            if box[i,j] != box[-i,-j]:
                box[i,j] = np.conjugate(box[-i,-j])
            else:
                box[i,j] = np.real(box[i,j])
            
    return box

            
    

## Counts the size of the power spectrum array produced by a n x n box


def count_size(n):
    x = fftshift(np.random.normal(size=(n,n)))
    ps = power_spec(x)
    ps_size = len(ps)
    return ps_size
    


