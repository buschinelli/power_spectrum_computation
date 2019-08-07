
import numpy as np
from numpy.fft import fft, fft2, fftfreq, ifft, ifft2, fftshift, fftn, ifftn, ifftshift
import matplotlib.pyplot as plt

def k_3d(box,dx):
    N = box.shape[0]
    k_comp = fftshift(fftfreq(N,d=dx))
    ans = np.zeros((N,N,N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                ans[i,j,k] = np.sqrt((k_comp[i]**2+k_comp[j]**2 + k_comp[k]**2))
    return ans


def k_2d(box,dx):
    N = box.shape[0]
    k_comp = fftshift(fftfreq(N,d=dx))
    ans = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            ans[i,j] = np.sqrt((k_comp[i]**2+k_comp[j]**2))
            
    return ans


def sigma_hat(box,dx):
    ibox = (dx)**3*fftshift(fftn(fftshift(box)))

    N = box.shape[0]    
    dk = 1/(N*dx)
    
    #ik,jk,kk = arr_kindex[0],arr_kindex[1],arr_kindex[2]
    sigma_hat = np.zeros((2*N,2*N,2*N))
    
    
    for ik1 in range(-N//2,N//2):
        for jk1 in range(-N//2,N//2):
            for kk1 in range(-N//2,N//2):
                for ik2 in range(-N//2,N//2):
                    for jk2 in range(-N//2,N//2):
                        for kk2 in range(-N//2,N//2):
                            ik,jk,kk = ik1 - ik2,jk1 - jk2, kk1 - kk2
                            sigma_hat[ik,jk,kk] += np.abs(ibox[ik1,jk1,kk1]*np.conjugate(ibox[ik2,jk2,kk2]))**2
            
                            
                
                
                
    sigma_hat *= (dk)**3/(2*np.pi)**3
    
    return fftshift(sigma_hat)



def trunk_sigma(sigma_hat,bins,dx):
    N = sigma_hat.shape[0]
    sigma_z = np.zeros((bins))
    count = np.zeros((bins))
    k_dist = k_3d(sigma_hat,dx)
    k_max = np.amax(k_dist)
    delta_k = k_max/bins
    k_plot = (delta_k/2)*np.arange(1,(2*bins+1),2) 
    
    for k in range(N):
        for i in range(N):
            for j in range(N):
                z = np.sqrt((i-N//2)**2+(j-N//2)**2+(k-N//2)**2)/(N*dx)
                loc_bin = int(z // delta_k)
                if z >= k_max:
                    loc_bin = bins-1
                
                sigma_z[loc_bin] += sigma_hat[i,j,k]
                count[loc_bin] += 1
                
    return sigma_z/count, k_plot
                
        
                
