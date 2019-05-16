import numpy as np
import math as m
from numpy.fft import fft, fft2, fftfreq, ifft, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt



## Take a NxN box in real space, compute its FFT and take its absolute
## value.

def abs_field(box):
    fft_box = fftshift(fft2((box)))
    psi_sqr = np.abs(fft_box)

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






        
    
