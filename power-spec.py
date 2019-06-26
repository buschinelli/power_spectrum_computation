import numpy as np
import math as m
import cmath as cm
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt
import random as rand




## Takes an NxN box in real space, compute its FFT and take its absolute
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

## Locates in what entries you have complex numbers, and where would you find its conjugate. Elements output 0 when 
## there is no complex conjugate for those elements


def location_cc(complex_box):
    real_box = np.real(complex_box)
    N = int(np.sqrt(np.size(real_box)))

    zeros = np.zeros((N,N))

    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                for h in range(0,N):
                    if real_box[i,j] == real_box[k,h]:
                        if i != k or j != h:
                            zeros[i,j] = real_box[i,j]
                            zeros[k,h] = real_box[k,h]
    return zeros

## This is the main code: it uses as input a given power spectrum and compute some box that would have as power spec
## the array ps input that way

def find_box(ps):
    size = np.size(ps)
    N = count_dim(size)
    box = np.zeros((N,N),dtype='complex')
    sample_box = np.random.normal(size=(N,N))  
    ft_sample_box = fftshift(fft2(sample_box))
    loc = location_cc(ft_sample_box)
    size = 0
    K = fftshift(k_matrix(box))
    k_value = np.unique(K)
    for i in range(0,N):
        for j in range(0,N):
            if box[i,j] == loc[i,j]:
                for p in range(0,len(k_value)):    
                    number = np.random.normal(0,np.sqrt(ps[p]),size=(1,1))
                    box[i,j] = number[0]
                    
            else:
                for p in range(0,len(k_value)):
                    size = size + 1
                number_c = rand_complex(np.sqrt(ps[p])*size,2)
                box[i,j] = number_c[0]
                box[-i,-j] = number_c[1]
                size = 0
        

    


    return ifft2(box)
    
    
    


## Generates a 1 x elements matrix with random complex number followed by their complex conjugate. If elements is odd
## then the last complex number doesn't have its pair
        

def rand_complex(value,elements):
    ps = value/(2*elements)
    i = cm.sqrt(-1)
    z = np.zeros(elements,dtype='complex')
    cc = 0
    while cc < elements:
        a = np.random.normal(0,ps,size=(1,1))
        b = np.random.normal(0,ps,size=(1,1))
        z[cc] = a[0] + i*b[0]
        cc = cc + 1
        if cc == elements:
            return z
        z[cc] = a[0] - i*b[0]
        cc = cc + 1
        
            
    
    return z

## Counts the size of the power spectrum array produced by a n x n box


def count_size(n):
    x = fftshift(np.random.normal(size=(n,n)))
    ps = power_spec(x)
    ps_size = len(ps)
    return ps_size
    





