from functools import reduce
import numpy as np
from constants import *
from numpy import linalg as la
from green_function import *


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


def onsite(t, mu):
    return (2 * t - mu) * s_z


def hopping(t, delta):
    return -t * s_z +  1j * delta * s_y


kitaev_chain = Hamiltonian()
infinite_nanowire = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
infinite_nanowire[lat(0)] = onsite
infinite_nanowire[kwant.HoppingKind((1,), lat)] = hopping
finite_nanowire = kwant.Builder()
finite_nanowire.fill(infinite_nanowire, (lambda site: 0 <= site.pos[0] < 20), (0,))
infinite_nanowire = kwant.wraparound.wraparound(infinite_nanowire).finalized()
finite_nanowire = finite_nanowire.finalized()

J = 2.0
p1 = dict(t=J / 2, mu=-1 * J, B=J, delta=2 * J, alpha=J)
p2 = dict(t=J / 2, mu=-3 * J, B=J, delta=2 * J, alpha=J)

H1 = finite_nanowire.hamiltonian_submatrix(params=p1)
H2 = finite_nanowire.hamiltonian_submatrix(params=p2)

h1_k = lambda k_x: infinite_nanowire.hamiltonian_submatrix(params=dict(**p1, k_x=k_x))
h2_k = lambda k_x: infinite_nanowire.hamiltonian_submatrix(params=dict(**p2, k_x=k_x))

periods = np.linspace(0.2 / J, 1.6 / J, 100)
momenta = np.linspace(-np.pi, np.pi)

energies = calculate_finite_spectrum(periods, [H1, H2])
spectrum = np.array([calculate_bands(momenta, [h1_k, h2_k], T) for T in periods])


def plot(n):
    T = J * periods[n]

    plot_1 = holoviews.Path(
        (J * periods, energies),
        kdims=[r"Driving period $(JT)$", r"Quasi-energy $(ET)$"],
        label="Finite system",
    ).opts(plot={"xticks": 5, "yticks": pi_ticks})

    VLine = holoviews.VLine(T).opts(style={"color": "b", "linestyle": "--"})

    plot_2 = holoviews.Path(
        (momenta, spectrum[n]), kdims=["$k$", "$E_kT$"], label="Floquet bands"
    ).opts(plot={"xticks": pi_ticks, "yticks": pi_ticks, "aspect": "equal"})
    return plot_1 * VLine + plot_2


holoviews.HoloMap({n: plot(n) for n in np.arange(0, 100, 10)}, kdims=["n"]).collate()