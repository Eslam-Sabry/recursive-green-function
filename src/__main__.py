from green_function import *


var_mu = np.linspace(-5*t,5*t,10001)
mu = .0*t
F_img = np.zeros([len(var_mu)])
G_img = np.zeros([len(var_mu)])
for i in range(len(var_mu)):
    #u = - var_mu[i] * s_z
    #u = - mu * s_z
    #v = -t * s_z + 1j * delta * s_y
    
    #u =  (-mu-2*t) * np.kron(s_z, s_0) + delta * np.kron(s_y, s_y) + B * np.kron(s_z, s_x)
    #v = -t * np.kron(s_z, s_0) - 1j * a * np.kron(s_0, s_z)

    u =  -omega * np.kron(j_z,np.kron(s_0, s_0)) + (-mu-2*t) * np.kron(j_0,np.kron(s_z, s_0)) + delta * np.kron(j_0, np.kron(s_y, s_y)) + B * np.kron(j_x, np.kron(s_z, s_y-s_x))
    v = -t * np.kron(j_0, np.kron(s_z, s_0)) - 1j * a * np.kron(j_0, np.kron(s_0, s_z))

    
    ##u = (2 * t - mu) * np.kron(s_z,s_0) + B * np.kron(s_0,s_z) + delta * np.kron(s_x,s_0)
    ##v = -t * np.kron(s_z,s_0) + 0.5 * 1j * a * np.kron(s_z,s_x)

    h = Hamiltonian(12, N, u, v)
    #G_img[i] = -1*np.array([left_green_function(0.00001j*t,h)[0,0].imag])
    #F_img[i] = -1*np.array([left_green_function(0.00001j*t,h)[0,1].imag])
    G_img[i] = np.trace(left_green_function(var_mu[i]-0.00001j*t,h).imag)
plt.title("Green function at the edge of semi-infinite Kitaev chain")
plt.plot(var_mu/t,G_img, label='G_img')
#plt.plot(var_mu/t,F_img,label='F_img')
#plt.plot(var_mu/t,G_img-F_img,label = 'G_img-F_img')
plt.ylabel('Green function')
plt.xlabel('$\mu/t$')
#plt.yscale('log')
plt.legend()
plt.show()

print(np.trace(left_green_function(-0.00001j*t,h).imag))
print(np.trace(left_green_function(t-0.00001j*t,h).imag))

