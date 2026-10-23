"""Check the results of the bifurcation example.

Run from this folder after the solver has finished:

    python3 check_results.py

Each output file <segment>_<field>.dat holds one row per finite element node
(first row = segment inlet, last row = segment outlet) and one column per
saved time step.
"""
import sys

import numpy as np

q0 = np.loadtxt("seg0_flow.dat")      # trachea
q1 = np.loadtxt("seg1_flow.dat")      # left bronchus
q2 = np.loadtxt("seg2_flow.dat")      # right bronchus
p1 = np.loadtxt("seg1_pressure.dat")

peak_flow = q0[-1].max()
junction_error = np.abs(q0[0] - q1[-1] - q2[-1]).max()
peak_pressure = np.abs(p1[0]).max()

print(f"Peak tracheal flow          : {peak_flow:.3f} mL/s   (expected ~392.7)")
print(f"Flow balance error at joint : {junction_error:.2e} mL/s (expected < 1e-6)")
print(f"Peak distal pressure        : {peak_pressure:.3f} dyn/cm2 (expected ~7.37)")

ok = (abs(peak_flow - 392.7) < 0.5
      and junction_error < 1e-6
      and abs(peak_pressure - 7.37) < 0.05)
print("OK" if ok else "UNEXPECTED RESULTS")
sys.exit(0 if ok else 1)
