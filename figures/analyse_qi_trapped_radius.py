"""Longer record of the trapped orbits of QINERTIA.md section 4: true max radius over time to t = 5000 (all orbits keep integrating after exit,
so only the radius of those that never crossed r = 5 is meaningful).  Usage: python3 figures/analyse_qi_trapped_radius.py"""
import sys, numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import quantised_inertia as Q
from lakes_of_wada.scattering import harmonic_section, integrate_exits
sysm = Q.hh_system(3, 5.0, "mond", 0.2)
y0, valid = harmonic_section(3, 80, 0.25); s0 = y0[valid]
lab, _, _ = integrate_exits(sysm, s0, 0.05, 200.0)
s = s0[lab == 0][:400].copy()
for k in range(int(5000/dt)):
    k1=f(s); k2=f(s+.5*dt*k1); k3=f(s+.5*dt*k2); k4=f(s+dt*k3); s = s+dt/6*(k1+2*k2+2*k3+k4)
    r = np.hypot(s[:,0], s[:,1]); rmax = np.maximum(rmax, r); exited |= r > 5.0
print("orbits", len(s), "t=5000 dt=0.02; exited(r>5) %.3f" % exited.mean())
print("true max r over time: percentiles 50/90/99/100", np.percentile(rmax[~exited],[50,90,99,100]).round(3))
print("saddle radius 1, saddle energy 1/6=0.1667")
