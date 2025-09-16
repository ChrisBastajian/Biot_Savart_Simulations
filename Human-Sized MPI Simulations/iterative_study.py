import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad, dblquad
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = 'browser'
#curve path of current:
R = 10 * 1e-2 # the radius of the coil, assuming 90th percentile to be no bigger than 20 (overestimation)
mu0 = 4 * np.pi * 1e-7
I = 10 #current intensity in A
h = 20 * 1e-2 #~20cm, assuming 90th percentile to be no bigger than 18 (overestimation

#finding thickness of litz wire:
data_points_coil = 100000
num_layers = 4 #4 layers of coils
coil_thickness = 3 * 1e-3

#getting number of turns:
num_turns_per_layer = h/coil_thickness
num_turns = num_turns_per_layer * num_layers
print(f"Number of turns: {num_turns}")

#Function that returns the integrand based on dimensions:
def get_integrand(radius, n_turns, height):
    phi, x, y, z = sp.symbols('phi,x,y,z')

    #recreate function of current path but in sympy to integrate:
    l = radius * sp.Matrix([sp.cos(phi), sp.sin(phi), (height/(radius * n_turns * 2 * np.pi))*phi])

    #distance away from path:
    r = sp.Matrix([x, y, z])

    #B(r) = ((mu0 * I)/(4piR)) * integral_along_curve\loop([(dl/dt x (r-l))/norm(r-l)^3)]dt
    #so we can simplify by defining r-l (vector of separation):
    difference = r-l
    #print(difference)

    #function to integrate (commented above):
    integrand = (mu0 * I / (4 * np.pi)) * sp.diff(l, phi).cross(difference)/ difference.norm()**3 #this will hold dB/dt for all 3 dimensions

    # need to lambdify all 3 functions into actual array representations of them:
    dBx = sp.lambdify([phi, x, y, z], integrand[0])
    dBy = sp.lambdify([phi, x, y, z], integrand[1])
    dBz = sp.lambdify([phi, x, y, z], integrand[2])

    return integrand, dBx, dBy, dBz

#function to get B from the integrand (integrate):
def B(xo, yo, zo, dBx, dBy, dBz, n_turns):
    return(np.array([
        quad(dBx, 0, n_turns * 2 * np.pi, args=(xo, yo, zo), limit=10000)[0],
        quad(dBy, 0, n_turns * 2 * np.pi, args=(xo, yo, zo), limit=10000)[0],
        quad(dBz, 0, n_turns * 2 * np.pi, args=(xo, yo, zo), limit=10000)[0],
    ]))

def get_inductance(radius, n_turns, height):
    mu_o = 4 * np.pi * 1e-7
    L = (np.pi * pow(n_turns,2) * pow(radius,2) * mu_o)/height
    return L

integrands, dBx, dBy, dBz = get_integrand(R, num_turns, h)
print(dBz(2*np.pi, 1, 1,1 )) #prints the z component at the given location
#calculate B at the center (0, 0, 0)
B_center = B(0, 0, h/2, dBx, dBy, dBz, num_turns)
print(B_center)


print(get_inductance(0.0175, 196, 0.07))
#Start of Iterative Process:
