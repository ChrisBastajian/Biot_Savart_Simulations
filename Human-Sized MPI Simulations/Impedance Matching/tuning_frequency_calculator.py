"""
This script will calculate the frequency of tuning for various capacitor values based on the agreed upon
coil design for the Tx coil for the test study
"""
import numpy as np

L= 0.35 *1e-3 #approx. 0.35mH

def get_tuning_frequency(C, L):
    #C = 1/(pow(w,2)*L)
    w = 1/np.sqrt(C*L)
    f = w/(2*np.pi)
    return f

Ceq_min = 45 * 1e-9
Ceq_max = 2.9 * 1e-6

print(f"Minimum Frequency: {get_tuning_frequency(Ceq_max, L)}",
      f"Maximum Frequency: {get_tuning_frequency(Ceq_min, L)}")

C_vals = np.array([0.125,0.25,0.5,1,1.5,2,3]) *1e-6
for C in C_vals:
    f = get_tuning_frequency(C, L)
    print(f"C = {C}, F = {f}")