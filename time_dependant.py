# This is to simulate an Ac field and see the change over time for one circular loop

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.vector import cross
from scipy.integrate import quad

t = np.linspace(0, 0.01, 100) #from 0 to 10 seconds
frequency = 1000 #Hz
w = 2*np.pi * frequency

I = 5 * np.cos(w*t) #5A

#plt.plot(t, I)
#plt.show()

#Loop paramaters:
R = 0.005 #50 mm
datapoints = 1000
mu0 = 4 * np.pi * 1e-7

#the radial path that current will follow:
phi = np.linspace(0, 2 * np.pi, datapoints)

l = R * np.array([np.cos(phi), np.sin(phi), np.zeros(len(phi))])
lx, ly, lz = l[0], l[1], l[2]

############################################################################
#using sympy to calculate the integral:
t, phi, x, y, z = sp.symbols('t, phi, x, y, z')

#recreate function of current path but in sympy to integrate:
l = R * sp.Matrix([sp.cos(phi), sp.sin(phi), 0])

#distance away from path:
r = sp.Matrix([x, y, z])

#B(r) = ((mu0 * I)/(4piR)) * integral_along_curve\loop([(dl/dt x (r-l))/norm(r-l)^3)]dt
#so we can simplify by defining r-l (vector of separation):
difference = r-l
#print(difference)

i = 5 * sp.cos(w*t)

#function to integrate (commented above):
integrand = (mu0 * i / (4 * np.pi)) * sp.diff(l, phi).cross(difference)/ difference.norm()**3 #this will hold dB/dt for all 3 dimensions

print(integrand)

#need to lambdify all 3 functions into actual array representations of them:
dBx = sp.lambdify([t, phi,x, y, z], integrand[0])
dBy = sp.lambdify([t, phi,x, y, z], integrand[1])
dBz = sp.lambdify([t, phi,x, y, z], integrand[2])

print(dBz(0, 2*np.pi, 1, 1,1 )) #prints the z component at the given location

#function to get B from the integrand (integrate):
def B(t, x, y, z):
    return(np.array([
        quad(dBx, 0, 2 * np.pi, args=(t, x, y, z), limit=10000)[0],
        quad(dBy, 0, 2 * np.pi, args=(t, x, y, z), limit=10000)[0],
        quad(dBz, 0, 2 * np.pi, args=(t, x, y, z), limit=10000)[0],
    ]))

H = np.linspace(0, 1, 100)
for j in range(100):
    n = j/10000
    H[j] = B(n, 0, 0, 0)[2]

plt.plot(H)
plt.show()