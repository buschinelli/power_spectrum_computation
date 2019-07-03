import numpy as np
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt
import random as rand

def power_spec(box,dx,bins):
    N = box.shape[0]
    ibox = (dx**2)*fftshift(fft2(fftshift(box)))
    ibox_abs_sq = np.abs(ibox)**2
    k_comp = fftshift(fftfreq(N,d=dx))
    K = 2*np.pi*np.sqrt(np.outer(np.ones(N),k_comp)**2 + np.outer(k_comp,np.ones(N))**2)
    unique_K = np.unique(K,return_counts=True)
    max_k = np.amax(unique_K[0])
    ps = np.zeros((bins))
    step = max_k/bins # Delta-k

    
    n_elements = np.zeros((bins))
    for i in range(-N//2,N//2):
        for j in range(-N//2,N//2):
            z = 2*np.pi*np.sqrt(i**2+j**2)/(N*dx)
            bin_num = int(z // step)-1
            ps[bin_num] += ibox_abs_sq[i,j]
            n_elements[bin_num] += 1
    
    ps /= n_elements
    
    return ps/((N*dx)**2)  



def find_box(ps,dx,N):
    bins = ps.shape[0]
    a = np.random.normal(size=(N,N))
    b = np.random.normal(size=(N,N))
    box = a + 1j*b
    k_comp = fftshift(fftfreq(N,d=dx))
    shift_box = fftshift(box)
    K = 2*np.pi*np.sqrt(np.outer(np.ones(N),k_comp)**2 + np.outer(k_comp,np.ones(N))**2)
    unique_K = np.unique(K,return_counts=True) 
    max_k = np.amax(unique_K[0])
    step = max_k/bins
    for i in range(-N//2,N//2):
        for j in range(-N//2,N//2):
            z = np.sqrt(i**2+j**2)/(N*dx)
            bin_num = int(z // step)-1
            box[i,j] *= np.sqrt((N*dx)**2*ps[bin_num]/2)
            
    
    for i in range(0,N):
        for j in range(0,N):
            if box[i,j] != box[-i,-j]:
                box[i,j] = np.conjugate(box[-i,-j])
            else:
                box[i,j] = np.real(box[i,j])
                
    
    return ifft2(box)/(dx**2)




## Tests

dx = 0.23
N = 200
bins = 80

box = np.random.normal(size=(N,N))
average = np.zeros((bins))
universes = 100
all_ps = np.zeros((bins))

for i in range(universes):
    ps = power_spec(box,dx,bins)
    all_ps = np.vstack((all_ps,ps))
    box = find_box(ps,dx,N)   
    
averages = np.mean(all_ps,axis=0)
plt.plot(averages)
plt.show()




    
    
