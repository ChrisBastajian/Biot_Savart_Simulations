"""
The axis system is such that the datum is positioned at the bottom of the Tx coil. This script calculates the turn
ratios needed between Rx and Cx in order to cancel the signal out when they are closest/ farthest to create a range
of number of turns in Cx per turn in Rx for tuning.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import sympy as sp
import plotly.graph_objects as go

z = 0.14 #m offset of Tx coil in z-dir
b = 0.01 #m
sq = 0.02 #m
n =  0.003 #m
r = 0.015 #m
l_supp = 0.02 #m
H = 0.07 #m
Hdisk = 0.035 #m
l_nmr = 0.03 # m (from top of Tx **COIL** to center of Rx coil)

k = 0.5*(H-r)
k_prime = k-n

L = z+b-sq #height of Cx with its holder part
j = L-(k_prime+Hdisk+r+n)

d_rx_cx = 2*n +r #distance between Cx and Rx coils (between centers)

Rx_pos = H - l_nmr
lowest_Cx = -(z-sq - k_prime - Hdisk - j - 0.5*r)
highest_Cx = H - l_nmr - r - 2*n

print(f"Lowest Cx position in z: {lowest_Cx}\n"
      f"Highest Cx position in z: {highest_Cx}\n"
      f"Rx position in z: {Rx_pos}\n"
      f"Distance between Rx and Cx when closest: {d_rx_cx}\n"
      f"Distance between Rx and Cx when furthest: {abs(lowest_Cx) + H - l_nmr}\n")

#curve path of current:
R = 0.042 # the radius of the curve itself (not the distance for each vector) in m
mu0 = 4 * np.pi * 1e-7
I = 5 #current intensity in A
n_turns = 636
max_height = 0.133

#finding thickness of litz wire:
total_thickness = 0.008 #mm
data_points_coil = 100000
num_coils = 4
coil_thickness = total_thickness/num_coils

def coil_arrays(num_coils, coil_thickness):
    phi = {}

    for n in range(num_coils):
        n_string = str(n)
        phi[n_string] = np.linspace(0, n_turns *2 * np.pi, round(data_points_coil/num_coils))

    return phi

a = coil_arrays(num_coils, coil_thickness)
#print(a['3']/a['1']) #to check that all arrays are the same
#print(len(a['1'])) #and that they have the length corresponding to dividing the datapoints wanted by the number of coils
#height = np.linspace(0, max_height, 100000)
#phi = np.linspace(0, n_turns * 2* np.pi, 100000)

def l(phi, R, coil_thickness, coil_number):
    R = R + coil_thickness * coil_number
    lengths = R * np.array([np.cos(phi), np.sin(phi), (max_height / (R* (n_turns * 2 * np.pi))) * phi]) #to create a circle with increments (loops)
    return lengths

#store dictionaries of lx, ly, and lz:
lx, ly, lz = {}, {}, {}

#all 3 dimensions for the curves:
for m in range(num_coils):
    m_string = str(m)
    lx[m_string], ly[m_string], lz[m_string] = l(a[m_string], R, coil_thickness, m + 1)

#Plot to check that the curves are all what we want:
plt.plot(lx['0'],ly['0'])
plt.plot(lx['1'],ly['1'])
plt.plot(lx['2'],ly['2'])
plt.plot(lx['3'], ly['3'])
plt.show()

#solve for integrals using sympy:
phi, x, y, z = sp.symbols('phi,x,y,z')

#recreate function of current path but in sympy to integrate:
l = R * sp.Matrix([sp.cos(phi), sp.sin(phi), (max_height/(R * n_turns * 2 * np.pi))*phi])

#distance away from path:
r = sp.Matrix([x, y, z])

#B(r) = ((mu0 * I)/(4piR)) * integral_along_curve\loop([(dl/dt x (r-l))/norm(r-l)^3)]dt
#so we can simplify by defining r-l (vector of separation):
difference = r-l
#print(difference)

#function to integrate (commented above):
integrand = (mu0 * I / (4 * np.pi)) * sp.diff(l, phi).cross(difference)/ difference.norm()**3 #this will hold dB/dt for all 3 dimensions

print(integrand)

#need to lambdify all 3 functions into actual array representations of them:
dBx = sp.lambdify([phi,x, y, z], integrand[0])
dBy = sp.lambdify([phi,x, y, z], integrand[1])
dBz = sp.lambdify([phi,x, y, z], integrand[2])

print(dBz(2*np.pi, 1, 1,1 )) #prints the z component at the given location

#function to get B from the integrand (integrate):
def B(x, y, z):
    return(np.array([
        quad(dBx, 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0],
        quad(dBy, 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0],
        quad(dBz, 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0],
    ]))

#calculate B at the center (0, 0, 0)
B_Rx = B(0, 0, Rx_pos)
print("B_Rx:", B_Rx[2])
B_Cx_low = B(0,0,lowest_Cx) #the sample is 45 mm from the top
print(f"B_Cx_low: {B_Cx_low[2]}")
B_Cx_high = B(0,0,highest_Cx)
print(f"B_Cx_high: {B_Cx_high[2]}")

#B_Cx = n B_Rx when Cx is at high position:
n = B_Rx/B_Cx_high
print(f"Num_turns high: {n[2]}")

n = B_Rx/B_Cx_low
print(f"Num_turns low: {n[2]}")