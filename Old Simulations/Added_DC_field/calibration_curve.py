import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from matplotlib import pylab
from sympy.vector import cross
from scipy.integrate import quad
import plotly.graph_objects as go

#curve path of current:
r = 0.042 # the radius of the curve itself (not the distance for each vector) in m
mu0 = 4 * np.pi * 1e-7
I = 0 #current intensity in A
n_turns = 626
max_height = 0.133

#finding thickness of litz wire:
#total_thickness = 0.008 #in m (we use #20 AWG = 0.81 mm)
data_points_coil = 100000
num_coils = 6
n_turns = n_turns/num_coils #number of turns per coil
#coil_thickness = total_thickness/num_coils
coil_thickness = 0.81 * 1e-3

def coil_arrays(num_coils):
    phi = {}

    for n in range(num_coils):
        n_string = str(n)
        phi[n_string] = np.linspace(0, n_turns *2 * np.pi, round(data_points_coil/num_coils))

    return phi

a = coil_arrays(num_coils)
#print(a['3']/a['1']) #to check that all arrays are the same
#print(len(a['1'])) #and that they have the length corresponding to dividing the datapoints wanted by the number of coils
#height = np.linspace(0, max_height, 100000)
#phi = np.linspace(0, n_turns * 2* np.pi, 100000)

def l(phi, r, coil_thickness, coil_number):
    r = r + coil_thickness * coil_number
    lengths = r * np.array([np.cos(phi), np.sin(phi), (max_height / (r* (n_turns * 2 * np.pi))) * phi]) #to create a circle with increments (loops)
    return lengths

#store dictionaries of lx, ly, and lz:
lx, ly, lz = {}, {}, {}

#all 3 dimensions for the curves:
for m in range(num_coils):
    m_string = str(m)
    lx[m_string], ly[m_string], lz[m_string] = l(a[m_string], r, coil_thickness, m + 1)

#Plot to check that the curves are all what we want:

try:
    fig1 = plt.figure(figsize=(16, 8))

    plt.subplot(1, 2, 1)
    plt.plot(lx['0'], ly['0'])
    plt.plot(lx['1'], ly['1'])
    plt.plot(lx['2'], ly['2'])
    plt.plot(lx['3'], ly['3'])
    plt.title("Current Path X_Y")
    plt.xlabel("X")
    plt.ylabel("Y")

    ax = fig1.add_subplot(122, projection='3d')

    ax.plot(lx['0'], ly['0'], lz['0'], label='Coil 1', color='red')
    ax.plot(lx['1'], ly['1'], lz['1'], label='Coil 2', color='black')
    ax.plot(lx['2'], ly['2'], lz['2'], label='Coil 3', color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.show()

except:
    print('not enough coils for that')



##########################################################################################################################

#solve for integrals using sympy:
phi, x, y, z = sp.symbols('phi,x,y,z')

#set libraries:
l, R, difference, integrands = {}, {}, {}, {}
dBx, dBy, dBz = {}, {}, {}


#distance from path:
d = sp.Matrix([x, y, z])

datapoints = 50
#set matrices:
I_array, B_SPIO_z = np.zeros(datapoints), np.zeros(datapoints)

for s in range(datapoints): #50 datapoints
    I = 5 * (s/datapoints) #current values up to 5A

    #recreate function of current path but in sympy to integrate:
    for k in range(num_coils):
        k_string = str(k)
        R[k_string] = r + (k+1) * coil_thickness
        l[k_string] = R[k_string] * sp.Matrix([sp.cos(phi), sp.sin(phi), (max_height/(R[k_string] * n_turns * 2 * np.pi))*phi])
        difference[k_string] = d - l[k_string]
        integrands[k_string] = (mu0 * I / (4 * np.pi)) * sp.diff(l[k_string], phi).cross(
            difference[k_string]) / difference[k_string].norm() ** 3  # this will hold dB/dt for all 3 dimensions
        # need to lambdify all 3 functions into actual array representations of them:
        dBx[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][0])
        dBy[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][1])
        dBz[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][2])


    #function to get B from the integrand (integrate):
    def B(x, y, z):
        Bx, By, Bz = 0, 0, 0
        for j in range(num_coils):
            j_string = str(j)
            Bx += quad(dBx[j_string], 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0]
            By += quad(dBy[j_string], 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0]
            Bz += quad(dBz[j_string], 0, n_turns * 2 * np.pi, args=(x, y, z), limit=10000)[0]

        return(np.array([Bx, By, Bz]))

    #calculate B at the center (0, 0, 0)
    #B_center = B(0, 0, max_height/2)
    #print("Magnetic field at the center:", B_center)
    B_SPIO = B(0,0, max_height - 0.045) #the sample is 45 mm from the top
    print("Magnetic field at location of SPIO (Z_component):", I, "    ", B_SPIO[2])

    B_SPIO_z[s] = B_SPIO[2]
    I_array[s] = I

figure = plt.figure()
ax = figure.add_subplot(111)

B_SPIO_z = B_SPIO_z * 1e+3 - 0.042 #turn into mT and account for Earth Magnetic field

ax.plot(I_array, B_SPIO_z, marker='o')
ax.set_title("Calibration Curve MPS", fontsize = "32")
ax.set_xlabel("I_DC [A]", fontsize = "18")
ax.set_ylabel("Magnetic field [mT]", fontsize="18")

#Set ticks to display:
x_ticks = np.linspace(0, 5,10)
y_ticks = np.linspace(0, max(B_SPIO_z), 10)

ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)

ax.grid()

z = np.polyfit(I_array, B_SPIO_z, 1)
p = np.poly1d(z)
pylab.plot(I_array,p(I_array),"r--")
# the line equation:
print ("y=%.6fx(%.6f)"%(z[0],z[1]))


plt.show()

