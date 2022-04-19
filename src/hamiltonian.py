import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import math
from constants import *

class Hamiltonian():
    def __init__(self,d,N,onsite,hopping):
        #d = degrees of freedom
        self.d = d
        #number of lattice points
        self.N = N
        #same site interaction
        self.onsite = onsite
        #nearest neighbors interaction
        self.hopping = hopping
    
    def lattice_hamiltonian(self,**kwarg):
        d, N, u, v = self.d, self.N, self.onsite(**kwarg), self.hopping(**kwarg)
        H=np.zeros([d*N,d*N],complex)
        for i in range(N-1):
            H[i*d:(i+1)*d,i*d:(i+1)*d] = u
            H[i*d:(i+1)*d,(i+1)*d:(i+2)*d] = v
            H[(i+1)*d:(i+2)*d,i*d:(i+1)*d] = np.conj(v.transpose())
        H[(N-1) * d:N * d, (N-1) * d:N * d] = u
        return H
    
    def k_space_hamiltonian(self,k,**kwarg):
        d, N, u, v = self.d, self.N, self.onsite(**kwarg), self.hopping(**kwarg)
        H_k = u + math.cos(k) * v.real + math.sin(k) * v.imag
        return H_k
    
    def plot_spectrum(self,**kwarg):
        var_k = np.linspace(-1*np.pi,1*np.pi,50)
        spectrum = []
        for i in range(len(var_k)):
            eval,evec=la.eigh([self.k_space_hamiltonian(var_k[i],**kwarg)])
            spectrum.append(eval[0])
        plt.title("Energy spectrum")
        plt.plot(var_k,spectrum, label='spectrum')
        #plt.plot(var_mu/t,G_img-F_img,label = 'G_img-F_img')
        plt.ylabel('E')
        plt.xlabel('$k$')
        plt.legend()
        plt.show()






def onsite(mu=mu,**kwarg):
    u = - mu * s_z
    return u


def hopping(t=t, delta=delta,**kwarg):
    v = -t * s_z + 1j * delta * s_y
    return  v

kitaev = Hamiltonian(-d,N,onsite,hopping)
#print(kitaev.k_space_hamiltonian(np.pi/2))
x = np.pi
p = dict(t=x, mu=-1*x,delta=0.3*x)
kitaev.plot_spectrum(**p)