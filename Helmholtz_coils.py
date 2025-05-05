#Note: in this code we consider two main coils, each consisting of their own multiple loops and ticknesses
# soyou can have 5 coils per main coil (that are not separated/ laying wire on top of itself to get the number of turns you want
# for each main coil

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.vector import cross
from scipy.integrate import quad
import plotly.graph_objects as go

#curve representing current path:
R = 0.07 #70 mm
mu0 = 4 * np.pi * 1e-7
I = 10 #current intensity in A

wire_thickness = 0.0012 #1.2mm

n_turns_per_coil = 11 #turns in z direction for each wire in x and y, itbasically shows the height
n_loops_xy_per_coil = 3 #there is 3 circular paths per coil each containing 11 turns going in z direction

height_single_loop = n_turns_per_coil * wire_thickness

data_points_total = 10000

def coil_arrays(num_coils): #gives the trajectory along which the current flows in the "phi" dimension and splits data points for each coil
    phi = {}

    for n in range(num_coils): #num_coils is the number of sub coils for each of the two main coils
        n_string = str(n)
        phi[n_string] = np.linspace(0, n_turns_per_coil *2 * np.pi, round(data_points_total/num_coils))

    return phi

a = coil_arrays(n_loops_xy_per_coil) #num_coils per main coil...
#print(a['1']/a['0'])

def l(phi, r, wire_thickness, coil_number, offset):
    r = r + wire_thickness *coil_number #radius of each loop that consists one of the two big coils

    lengths = r * np.array([np.cos(phi), np.sin(phi), offset/r + phi* (height_single_loop/(r * 2*np.pi *n_turns_per_coil)) ])
    return lengths


lx, ly, lz = {}, {}, {}
lower_offset = 0 #for first main coil (lower)
for i in range(n_loops_xy_per_coil):
    i_string = str(i)
    lx[i_string], ly[i_string], lz[i_string] = l(a[i_string], R, wire_thickness, i+1, lower_offset)

coils_offset = R #Since Helmholtz coils are separated by their equal radii
lx1, ly1, lz1 = {}, {}, {}
for j in range(n_loops_xy_per_coil): #trajectory for second main coil (upper)
    j_string = str(j)
    lx1[j_string], ly1[j_string], lz1[j_string] = l(a[j_string], R, wire_thickness, j+1, coils_offset)



#Plot to check that the curves are all what we want:

try:
    fig1 = plt.figure(figsize=(16, 8))

    plt.subplot(1, 2, 1)
    plt.plot(lx['0'], ly['0'])
    plt.plot(lx['1'], ly['1'])
    plt.plot(lx['2'], ly['2'])

    plt.plot(lx1['0'], ly1['0'])
    plt.plot(lx1['1'], ly1['1'])
    plt.plot(lx1['2'], ly1['2'])
    plt.title("Current Path X_Y")
    plt.xlabel("X")
    plt.ylabel("Y")

    ax = fig1.add_subplot(122, projection='3d')

    ax.plot(lx['0'], ly['0'], lz['0'], label='Coil 1', color='red')
    ax.plot(lx['1'], ly['1'], lz['1'], label='Coil 2', color='black')
    ax.plot(lx['2'], ly['2'], lz['2'], label='Coil 3', color='blue')

    ax.plot(lx1['0'], ly1['0'], lz1['0'], label='Coil 1', color='red')
    ax.plot(lx1['1'], ly1['1'], lz1['1'], label='Coil 2', color='black')
    ax.plot(lx1['2'], ly1['2'], lz1['2'], label='Coil 3', color='blue')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.show()

except:
    print('not enough coils for that')

##################################################################################################################
#Solve the Biot_Savart integral using sympy:
phi, x, y, z = sp.symbols('phi,x,y,z')

#set libraries:
l, r, difference, integrands = {}, {}, {}, {}
dBx, dBy, dBz = {}, {}, {}
dBx1, dBy1, dBz1 = {}, {}, {}

#distance from path:
d = sp.Matrix([x, y, z])

offset = {}
offset['0'], offset['1'] = 0, R #first coil has no offset, but second coil has an offset of R

# current path but in sympy:
for k in range(n_loops_xy_per_coil):
    k_string = str(k)
    r[k_string] = R + (k+1) * wire_thickness
    l[k_string] = r[k_string] * sp.Matrix([sp.cos(phi), sp.sin(phi), offset['0']/r[k_string] +
                                           (height_single_loop / (r[k_string] * n_turns_per_coil * 2 * np.pi)) * phi])
    difference[k_string] = d - l[k_string]
    integrands[k_string] = (mu0 * I / (4 * np.pi)) * sp.diff(l[k_string], phi).cross(
        difference[k_string]) / difference[k_string].norm() ** 3  # this will hold dB/dt for all 3 dimensions
    dBx[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][0])
    dBy[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][1])
    dBz[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][2])

for k in range(n_loops_xy_per_coil):
    k_string = str(k)
    r[k_string] = R + (k + 1) * wire_thickness
    l[k_string] = r[k_string] * sp.Matrix([sp.cos(phi), sp.sin(phi), offset['1'] / r[k_string] +
                                           (height_single_loop / (
                                                       r[k_string] * n_turns_per_coil * 2 * np.pi)) * phi])
    difference[k_string] = d - l[k_string]
    integrands[k_string] = (mu0 * I / (4 * np.pi)) * sp.diff(l[k_string], phi).cross(
        difference[k_string]) / difference[
                               k_string].norm() ** 3  # this will hold dB/dt for all 3 dimensions
    dBx1[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][0])
    dBy1[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][1])
    dBz1[k_string] = sp.lambdify([phi, x, y, z], integrands[k_string][2])


def B(x, y, z):
    Bx0, By0, Bz0 = 0, 0, 0
    Bx1, By1, Bz1 = 0, 0, 0
    Bx, By, Bz = 0, 0, 0
    for j in range(n_loops_xy_per_coil):
        j_string = str(j)
        Bx0 += quad(dBx[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]
        By0 += quad(dBy[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]
        Bz0 += quad(dBz[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]
        Bx1 += quad(dBx1[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]
        By1 += quad(dBy1[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]
        Bz1 += quad(dBz1[j_string], 0, n_turns_per_coil * 2 * np.pi, args=(x, y, z), limit=10000)[0]

        Bx = Bx0 + Bx1
        By = By0 + By1
        Bz = Bz0 + Bz1

    return(np.array([Bx, By, Bz]))

B_center = B(0, 0, R/2)
print(B_center)


#meshgrid to display:
i = np.linspace(-2*R, 2*R, 10)
j = np.linspace(0, offset['1'] + 2 * height_single_loop, 10)
xi,yi,zi = np.meshgrid(i,i,j)

#get the whole field by vectorizing the result of calling B with the meshgrid elements:
B_field = np.vectorize(B, signature='(),(),()->(n)')(xi,yi,zi) #the signature in this
                                    # case is about joining all three inputs into one vector source (n)
Bx = B_field[:,:,:,0]
By = B_field[:,:,:,1]
Bz = B_field[:,:,:,2]

#plot everything
data = go.Cone(x=xi.ravel(), y = yi.ravel(), z=zi.ravel(), u =Bx.ravel(),
               v=By.ravel(), w = Bz.ravel())
fig = go.Figure(data= data)
fig.add_scatter3d(x=lx['0'], y=ly['0'], z=lz['0'], mode='lines', line=dict(color= 'black'))
fig.add_scatter3d(x=lx['1'], y=ly['1'], z=lz['1'], mode='lines', line=dict(color= 'black'))
fig.add_scatter3d(x=lx['2'], y=ly['2'], z=lz['2'], mode='lines', line=dict(color= 'black'))
fig.add_scatter3d(x=lx1['0'], y=ly1['0'], z=lz1['0'], mode='lines', line=dict(color= 'black'))
fig.add_scatter3d(x=lx1['1'], y=ly1['1'], z=lz1['1'], mode='lines', line=dict(color= 'black'))
fig.add_scatter3d(x=lx1['2'], y=ly1['2'], z=lz1['2'], mode='lines', line=dict(color= 'black'))

fig.show()




# Plotting using matplotlib in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the current path
ax.plot(lx['0'], ly['0'], lz['0'], label='Coil 1', color='red')
ax.plot(lx['1'], ly['1'], lz['1'], label='Coil 2', color='black')
ax.plot(lx['2'], ly['2'], lz['2'], label='Coil 3', color='blue')
ax.plot(lx1['0'], ly1['0'], lz1['0'], label='Coil 1', color='red')
ax.plot(lx1['1'], ly1['1'], lz1['1'], label='Coil 2', color='black')
ax.plot(lx1['2'], ly1['2'], lz1['2'], label='Coil 3', color='blue')

# Plot the magnetic field vectors
ax.quiver(xi, yi, zi, Bx, By, Bz, length=0.01, normalize=True, color='blue')

# Set labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Magnetic Field Vectors and Current Path')

plt.legend()
plt.show()