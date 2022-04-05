import numpy as np


#t = 2
t = np.pi
delta = 0.50*t
mu = 0.0*t
B = 0.30*t
a = 1.0*t
omega = 0.1*t
A = 0.01*t
d = 2
N = 25

s_x = np.array([[0,1],[1,0]])
s_y = np.array([[0,-1j],[1j,0]])
s_z = np.array([[1,0],[0,-1]])
s_0 = np.eye(2)
j_x = np.array([[0,1,0],[1,0,1],[0,1,0]])
j_z = np.array([[1,0,0],[0,0,0],[0,0,-1]])
j_0 = np.eye(3)