import numpy as np
from scipy import integrate

# Note: t0 is required for the odeint function, though it's not used here.
def lorentz_deriv((x, y, z), t0, sigma=10., beta=8./3, rho=28.0):
	#Compute the time-derivative of a Lorenz system."""
	return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]

x0 = [1, 1, 1]  # starting vector
N_trajectories = 25

np.random.seed(1)
x0 = -15 +30 * np.random.random((N_trajectories,3))
t = np.linspace(0, 4, 1000)  # one thousand time steps
x_t = np.asarray([integrate.odeint(lorentz_deriv, x0i, t) for x0i in x0])
print(type(x_t))
print(np.shape(x_t))
x_t=np.asarray(x_t)
print(type(x_t))
