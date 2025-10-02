"""
This code calculates variations of capacitors in series and checks which ones fit within the tuning frequency
target range. The capacitors in stock are considered to be limited to 6 for each type (2 boxes) though this can
be changed as needed
"""
import numpy as np
import itertools

L= 0.35 *1e-3 #approx. 0.35mH
capacitors = [
    1000e-3, 4700e-3, 10e-3, 47e-3, 0.068e-3, 0.1, 0.22, 0.47, 0.68,
    1.0, 0.047, 0.1, 0.22, 0.33, 0.47, 0.22, 0.47, 0.68, 1.0,
    1.0, 2.2, 4.7, 6.8
] #uF

capacitors = [c*1e-6 for c in capacitors] #converting to F

stock = {c:6 for c in capacitors}  #max 6 for each capacitance value (if we order another box)
V_rating = 300# Vrms
target_range = (45e-9, 2.9e-6)  # [45n,2.9u]F
V_required = [800, 1500, 2000] #Vrms

def voltage_divider_series(C_values, V_total, freq):
    w = 2 * np.pi * freq
    impedances = [1/(w*c) for c in C_values]  #Xc for each
    Z_total = sum(impedances)
    voltages = [V_total * (Zc/Z_total) for Zc in impedances] #voltage divider for each
    return voltages

def get_tuning_frequency(C, L=0.35*1e-3):
    #C = 1/(pow(w,2)*L)
    w = 1/np.sqrt(C*L)
    f = w/(2*np.pi)
    return f

def get_series_capacitance(C_values, n=1, equal=False):
    Ceq_inv = sum(1/c for c in C_values)
    Ceq = 1/Ceq_inv if Ceq_inv >0 else 0
    return Ceq

def get_parallel_capacitance(C_values):
    return sum(C_values)

valid_solutions = []
global voltages
#Optimization loop:
for n_series in range(2, 8):
    #Generating all combinations of capacitor values
    for combo in itertools.combinations_with_replacement(capacitors, n_series):
        Ceq_series = get_series_capacitance(combo)

        #Computing stress voltages at chosen frequency
        freq_check = get_tuning_frequency(Ceq_series)
        if freq_check<10000:
            voltages = voltage_divider_series(combo, V_required[0], freq_check)
        elif 10000<freq_check<20000:
            voltages = voltage_divider_series(combo, V_required[1], freq_check)
        elif freq_check>20000:
            voltages = voltage_divider_series(combo, V_required[2], freq_check)
        if max(voltages) > V_rating:
            continue

        #Trying parallel equivalent (adding to current capacitance)
        max_strings = min(stock[c] // combo.count(c) for c in set(combo))
        for n_parallel in range(1, max_strings+1):
            Ceq_total = get_parallel_capacitance([Ceq_series]*n_parallel)

            if target_range[0] <= Ceq_total <= target_range[1]:
                f_res = get_tuning_frequency(Ceq_total)
                valid_solutions.append({
                    "combo": combo,
                    "N_series": n_series,
                    "N_parallel": n_parallel,
                    "Ceq": Ceq_total,
                    "f_resonance": f_res,
                    "voltages": voltages
                })

valid_solutions = sorted(valid_solutions, key=lambda x: x["Ceq"]) #sorting results by Ceq

for sol in valid_solutions:
    combo_str = " + ".join(f"{c*1e6:.3f}uF" for c in sol['combo'])
    volt_str = ", ".join(f"{v:.1f}V" for v in sol['voltages'])
    print(f"{sol['N_parallel']}x [{combo_str}] in series "
          f"=> Ceq={sol['Ceq']*1e9:.1f} nF, "
          f"f_res={sol['f_resonance']/1000:.1f} kHz, "
          f"V_caps=[{volt_str}]")