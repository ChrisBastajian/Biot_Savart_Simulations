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

C=0.186 *1e-6
f = get_tuning_frequency(C, L)
print(f"C = {C}, F = {f}")

def get_tuning_capacitance(f, Leq):
    w = 2*np.pi*f
    return 1/(pow(w,2)*L)

C = get_tuning_capacitance(20000, L)
print(f"Tuning Capacitance: {C}")