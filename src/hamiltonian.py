import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import math
#from constants import *

class Hamiltonian():
    def __init__(self,d,N,u,v):
        #d = degrees of freedom
        self.d = d
        #number of lattice points
        self.N = N
        #same site interaction
        self.u = u
        #nearest neighbors interaction
        self.v = v
    
    def lattice_hamiltonian(self):
        d, N, u, v = self.d, self.N, self.u, self.v
        H=np.zeros([d*N,d*N],complex)
        for i in range(N-1):
            H[i*d:(i+1)*d,i*d:(i+1)*d] = u
            H[i*d:(i+1)*d,(i+1)*d:(i+2)*d] = v
            H[(i+1)*d:(i+2)*d,i*d:(i+1)*d] = np.conj(v.transpose())
        H[(N-1) * d:N * d, (N-1) * d:N * d] = u
        return H
    
    def k_space_hamiltonian(self,k):
        d, N, u, v = self.d, self.N, self.u, self.v
        H_k = u + math.cos(k) * v.real + math.sin(k) * v.imag
        return H_k
    
    def plot_spectrum(self):
        var_k = np.linspace(-np.pi,np.pi,50)
        spectrum = []
        for i in range(len(var_k)):
            eval,evec=la.eigh([self.k_space_hamiltonian(k=var_k[i])])
            spectrum.append(eval[0])
        plt.title("Energy spectrum")
        plt.plot(var_k,spectrum, label='spectrum')
        #plt.plot(var_mu/t,G_img-F_img,label = 'G_img-F_img')
        plt.ylabel('E')
        plt.xlabel('$k$')
        plt.legend()
        plt.show()

