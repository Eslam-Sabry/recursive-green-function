from green_function import *


var_mu = np.linspace(-0*t,3*t,5000)
mu = 0.0*t
F_img = np.zeros([len(var_mu)])
G_img = np.zeros([len(var_mu)])
for i in range(len(var_mu)):
    u = - var_mu[i] * s_z
    #u = - mu * s_z
    v = -t * s_z + 1j * delta * s_y
    h = Hamiltonian(d, N, u, v)
    G_img[i] = -1*np.array([left_green_function(0.00001j*t,h)[0,0].imag])
    F_img[i] = -1*np.array([left_green_function(0.00001j*t,h)[0,1].imag])
    #G_img[i] = np.trace(left_green_function(var_mu[i]-0.00001j*t,h).imag)
plt.title("Green function at the edge of semi-infinite Kitaev chain")
plt.plot(var_mu/t,G_img, label='G_img')
plt.plot(var_mu/t,F_img,label='F_img')
#plt.plot(var_mu/t,G_img-F_img,label = 'G_img-F_img')
plt.ylabel('Green function')
plt.xlabel('$\mu/t$')
plt.yscale('log')
plt.legend()
plt.show()
