from functools import reduce
import numpy as np
from constants import *
from numpy import linalg as la
from green_function import *


class Floquet(Hamiltonian):
    def __init__(self, d, N, u, v, T,):
        super().__init__(d, N, u, v)