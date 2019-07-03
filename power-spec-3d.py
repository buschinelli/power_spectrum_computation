import numpy as np
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt
import random as rand

def power_spec_3d(box,dx,bins):
    N = box.shape[0]
    ibox = (dx**3)*fftshift(fftn(fftshift(box)))
    ibox_abs_sq = np.abs(ibox)**2
    k_comp = fftshift(fftfreq(N,d=dx))
    K = 2*np.pi*k_3d(box,dx)
    unique_K = np.unique(K,return_counts=True)
    max_k = np.amax(unique_K[0])
    ps = np.zeros((bins))
    step = max_k/bins # Delta-k

    
    n_elements = np.zeros((bins))
    for i in range(-N//2,N//2):
        for j in range(-N//2,N//2):
            for k in range(-N//2,N//2):
                z = 2*np.pi*np.sqrt((i**2+j**2+k**2))/(N*dx)
                bin_num = int(z // step) -1
                ps[bin_num] += ibox_abs_sq[i,j,k]
                n_elements[bin_num] += 1

    ps /= n_elements
    
    return ps/((N*dx)**3)  



def find_box_3d(ps,dx,N):
    bins = ps.shape[0]
    a = np.random.normal(size=(N,N,N))
    b = np.random.normal(size=(N,N,N))
    box = a + 1j*b
    k_comp = fftshift(fftfreq(N,d=dx))
    shift_box = fftshift(box)
    K = 2*np.pi*k_3d(box,dx)
    unique_K = np.unique(K,return_counts=True) 
    max_k = np.amax(unique_K[0])
    step = max_k/bins
    for i in range(-N//2,N//2):
        for j in range(-N//2,N//2):
            for k in range(-N//2,N//2):
                z = np.sqrt(i**2+j**2+k**2)/(N*dx)
                bin_num = int(z // step)
                box[i,j,k] *= np.sqrt((N*dx)**3*ps[bin_num]/2)
            
    
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                if box[i,j,k] != box[-i,-j,-k]:
                    box[i,j,k] = np.conjugate(box[-i,-j,-k])
                else:
                    box[i,j,k] = np.real(box[i,j,k])
                
    
    return ifftn(box)/(dx**3)


def k_3d(box,dx):
    N = box.shape[0]
    k_comp = fftshift(fftfreq(N,d=dx))
    ans = np.zeros((N,N,N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                ans[i,j,k] = np.sqrt((k_comp[i]**2+k_comp[j]**2 + k_comp[k]**2))
    return ans
    


dx = 0.23
N = 70
bins = 50

box = np.random.normal(size=(N,N,N))
average = np.zeros((bins))
universes = 15
all_ps = np.zeros((bins))

for i in range(universes):
    ps = power_spec_3d(box,dx,bins)
    all_ps = np.vstack((all_ps,ps))
    box = find_box_3d(ps,dx,N)   
    
averages = np.mean(all_ps,axis=0)
plt.plot(averages)
plt.show()


