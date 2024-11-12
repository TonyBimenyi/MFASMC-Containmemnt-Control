import numpy as np
import matplotlib.pyplot as plt

#Defining Parameters

beta = -0.2
rho = 1
eta = 1
lamda = 50
mu = 1
epsilon = 10**-5
alpha1 = 1
alpha2 = 1.5
T = 0.1
gamma1 = 0.15
gamma2 = 0.15
gamma3 = 0.45
gamma4 = 0.45
rT = 1024
m = 350
L = 100

# Generate constant desired trajectory
yd = 0.5 * (1 + np.sign(np.sin(np.linspace(0, 2 * np.pi, L + 1))))  # Square wave



