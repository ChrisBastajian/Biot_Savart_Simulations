import numpy as np

L= 0.35 *1e-3 #approx. 0.35mH
capacitors = [
    1000e-3, 4700e-3, 10e-3, 47e-3, 0.068e-3, 0.1, 0.22, 0.47, 0.68,
    1.0, 0.047, 0.1, 0.22, 0.33, 0.47, 0.22, 0.47, 0.68, 1.0,
    1.0, 2.2, 4.7, 6.8
] #uF

capacitors = [c*1e-6 for c in capacitors] #converting to F

stock = {c:3 for c in capacitors}  # max 3 for each capacitance value
V_rating = 300# Vrms
target_range = (45e-9, 180e-9)  # [45,180] nF
V_required = 2000 #Vrms

def get_tuning_frequency(C, L=0.35*1e-3):
    #C = 1/(pow(w,2)*L)
    w = 1/np.sqrt(C*L)
    f = w/(2*np.pi)
    return f

def get_series_capacitance(C=None, n=1, equal=False):
    global Ceq
    if C is None:
        C = []
    if equal: #if there are many capacitors of the same value
        Ceq = C/n
        return Ceq
    else:
        for i in range(len(C)):
            C_current = C[i]
            if i == 0:
                Ceq = C_current
            else:
                Ceq = (C_current * Ceq)/(C_current+Ceq)
        return Ceq

def get_parallel_capacitance(C1, C2, n):
    if C1==C2:
        Ceq = C1*n
    else:
        Ceq = C1+C2
    return Ceq

valid_solutions = []

# Loop over capacitor types
for c in capacitors:
    for n_series in range(2, stock[c]+1):  # need ≥7 for 2kVrms
        # Ceq for series string
        Ceq_series = get_series_capacitance(C=c, n=n_series, equal=True)
        V_total = V_rating * n_series
        if V_total < V_required:
            continue

        # How many such strings can we build?
        max_strings = stock[c] // n_series
        for n_parallel in range(1, max_strings+1):
            Ceq_total = get_parallel_capacitance([Ceq_series]*n_parallel)

            if target_range[0] <= Ceq_total <= target_range[1]:
                f_res = get_tuning_frequency(Ceq_total)
                valid_solutions.append({
                    "C_single": c,
                    "N_series": n_series,
                    "N_parallel": n_parallel,
                    "Ceq": Ceq_total,
                    "V_rating_total": V_total,
                    "f_resonance": f_res
                })

valid_solutions = sorted(valid_solutions, key=lambda x: x["Ceq"]) #sorting results by Ceq

for sol in valid_solutions:
    print(f"{sol['N_parallel']} string(s) of {sol['N_series']}x {sol['C_single']*1e6:.3f} uF caps "
          f"=> Ceq={sol['Ceq']*1e9:.1f} nF, V={sol['V_rating_total']} Vrms, "
          f"f_res={sol['f_resonance']/1000:.1f} kHz")