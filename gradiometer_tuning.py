"""
The axis system is such that the datum is positioned at the bottom of the Tx coil. This script calculates the turn
ratios needed between Rx and Cx in order to cancel the signal out when they are closest/ farthest to create a range
of number of turns in Cx per turn in Rx for tuning.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
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