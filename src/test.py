from functools import reduce
import numpy as np
from constants import *
from green_function import *
#import holoviews
from scipy import linalg as la
import matplotlib as mp
plt.style.use('seaborn')		# Setting the plotting style
mp.rcParams['figure.figsize'] = (15, 10)  # Setting the size of the plots


def evolution_operator(hamiltonians, T):
    n = len(hamiltonians)
    exps = [la.expm(-1j * h * T / n) for h in hamiltonians]
    return reduce(np.dot, exps)


def calculate_finite_spectrum(periods, hamiltonians):
    energies = []
    for T in periods:
        U = evolution_operator(hamiltonians, T)
        phases = np.angle(la.eigvals(U))
        phases = np.sort(np.abs(phases))
        ev = np.sort([(-1) ** n * val for n, val in enumerate(phases)])
        energies.append(ev)
    return np.array(energies).real


def calculate_bands(momenta, hamiltonians_k, T):
    energies = []
    for k in momenta:
        hamiltonians = [h_k(k) for h_k in hamiltonians_k]
        U = evolution_operator(hamiltonians, T)
        phases = np.angle(la.eigvals(U))
        phases = np.sort(np.abs(phases))
        ev = np.sort([(-1) ** n * val for n, val in enumerate(phases)])
        energies.append(ev)
    return np.array(energies).real

'''
def onsite(mu=mu):
    return -mu * s_z


def hopping(t=t, delta=delta, A=A):
    return -t * s_z +  1j * delta * s_y + 1j * A * s_z


#var_time = dict()
periods = np.linspace(0.2 / t, 1.6 / t, 100)
momenta = np.linspace(-np.pi, np.pi)
#omegas = 2 * np.pi / periods
#var_time = [np.linspace(0,T,11) for T in periods]

h_1 = Hamiltonian(2,100,onsite(),hopping(A=0))
h_2 = Hamiltonian(2,100,onsite(),hopping(A=0.5*t))


energies = calculate_finite_spectrum(periods, [h_1.lattice_hamiltonian(), h_2.lattice_hamiltonian()])
spectrum = np.array([calculate_bands(momenta, [h_1.k_space_hamiltonian, h_2.k_space_hamiltonian], T) for T in periods])

pi_ticks = [
    (-np.pi, r"$-\pi$"),
    (-np.pi / 2, r"$-\pi/2$"),
    (0, r"$0$"),
    (np.pi / 2, r"$\pi/2$"),
    (np.pi, r"$\pi$"),
]

def plot(n):
    T = t * periods[n]

    plot_1 = holoviews.Path(
        (t * periods, energies),
        kdims=[r"Driving period $(JT)$", r"Quasi-energy $(ET)$"],
        label="Finite system",
    ).opts(plot={"xticks": 5, "yticks": pi_ticks})

    VLine = holoviews.VLine(T).opts(style={"color": "b", "linestyle": "--"})

    plot_2 = holoviews.Path(
        (momenta, spectrum[n]), kdims=["$k$", "$E_kT$"], label="Floquet bands"
    ).opts(plot={"xticks": pi_ticks, "yticks": pi_ticks, "aspect": "equal"})
    return plot_1 * VLine + plot_2

holoviews.extension('matplotlib')


#holoviews.HoloMap({n: plot(n) for n in np.arange(0, 100, 10)}, kdims=["n"]).collate()

'''




'''
def onsite(t=t, mu=mu, B=B, delta=delta):
    return (2 * t - mu) * np.kron(s_z,s_0) + B * np.kron(s_0,s_z) + delta * np.kron(s_x,s_0)


def hopping(t=t, alpha=a):
    return -t * np.kron(s_z,s_0) + 0.5 * 1j * alpha * np.kron(s_z,s_x)

J = 2.0
p1 = {'t':J / 2, 'mu':-1 * J, 'B':J, 'delta':2 * J, 'alpha':J}
p2 = {'t':J / 2, 'mu':-3 * J, 'B':J, 'delta':2 * J, 'alpha':J}



periods = np.linspace(0.2 / J, 1.6 / J, 100)
momenta = np.linspace(-np.pi, np.pi)

h_1 = Hamiltonian(4,N,onsite(t = p1['t'],mu=p1['mu'],B=p1['B'],delta=p1['delta']),hopping(t=p1['t'],alpha=p1['alpha']))
h_2 = Hamiltonian(4,N,onsite(t = p2['t'],mu=p2['mu'],B=p2['B'],delta=p2['delta']),hopping(t=p2['t'],alpha=p2['alpha']))


energies = calculate_finite_spectrum(periods, [h_1.lattice_hamiltonian(), h_2.lattice_hamiltonian()])
spectrum = np.array([calculate_bands(momenta, [h_1.k_space_hamiltonian, h_2.k_space_hamiltonian], T) for T in periods])


plt.title("Energy Spectrum as a function of $T$ (N = 25)")
for i in range(2*N):
    plt.plot(periods,spectrum[:,i])
plt.ylabel('Energy')
plt.xlabel('$\mu/t$')
plt.show()


plt.title("Energy Spectrum as a function of $T$ (N = 25)")
for i in range(4*N):
    plt.plot(periods,energies[:,i])
plt.ylabel('Energy')
plt.xlabel('$\mu/t$')
plt.show()

'''




def onsite(mu=mu):
    return -mu * s_z


def hopping(t=t, delta=delta):
    return -t * s_z +  1j * delta * s_y


#var_time = dict()
#periods = np.linspace(0.2 / t, 1.6 / t, 100)
periods = np.linspace(0.2 / t, 3 / t, 200)
momenta = np.linspace(-np.pi, np.pi,100)
#omegas = 2 * np.pi / periods
#var_time = [np.linspace(0,T,11) for T in periods]
N = 50
h_1 = Hamiltonian(d,N,onsite(mu=1*t),hopping())
h_2 = Hamiltonian(d,N,onsite(mu=3*t),hopping())


energies = calculate_finite_spectrum(periods, [h_1.lattice_hamiltonian(), h_2.lattice_hamiltonian()])
#spectrum = np.array([calculate_bands(momenta, [h_1.k_space_hamiltonian, h_2.k_space_hamiltonian], T) for T in periods])

#T = 0.6

#interesting_bulk_spectrum = calculate_bands(momenta, [h_1.k_space_hamiltonian, h_2.k_space_hamiltonian],T= T)

"""
plt.title(f"Energy Spectrum as a function of the period $T$ (N = {N} sites)")
for i in range(200):
    plt.plot(periods,spectrum[:,i])
plt.ylabel('Energy')
plt.xlabel('$\mu/t$')
plt.show()
"""

plt.title(f"Energy Spectrum as a function of the period $T$ (N = {N} sites)")
for i in range(d*N):
    plt.plot(periods,energies[:,i]/np.pi)
plt.ylabel('Energy/pi')
plt.xlabel('$T$, periodic chemical potential case')
#plt.savefig('mu case-trivial.png')
plt.show()
'''
plt.title(f"Bulk system's Energy Spectrum as a function of momentum $k$ at T = {T}")
plt.plot(momenta,interesting_bulk_spectrum/np.pi)
plt.ylabel('Energy/pi')
plt.xlabel('$k$, periodic chemical potential case')
#plt.savefig('bulk spectrum.png')
plt.show()
'''








'''
def onsite(mu=mu):
    return -mu * s_z


def hopping(t=t, delta=delta,A=0):
    return ( -t * np.exp(1j*A)) * s_z +  1j * delta * s_y

A = 0.1*t
h_1 = Hamiltonian(d,N,onsite(mu=1*t),hopping(A = -A))
h_2 = Hamiltonian(d,N,onsite(mu=1*t),hopping(A = A))




energies = calculate_finite_spectrum(periods, [h_1.lattice_hamiltonian(), h_2.lattice_hamiltonian()])
#spectrum = np.array([calculate_bands(momenta, [h_1.k_space_hamiltonian, h_2.k_space_hamiltonian], T) for T in periods])
'''

'''
plt.title(f"Energy Spectrum as a function of the period $T$ (N = {N} sites)")
for i in range(200):
    plt.plot(periods,spectrum[:,i])
plt.ylabel('Energy')
plt.xlabel('$\mu/t$')
plt.show()

'''

'''
plt.title(f"Energy Spectrum as a function of the period $T$ (N = {N} sites)")
for i in range(d*N):
    plt.plot(periods,energies[:,i]/np.pi)
plt.ylabel('Energy/pi')
plt.xlabel('$T$, periodic external field case')
plt.savefig('field case-topological.png')
#plt.show()
'''