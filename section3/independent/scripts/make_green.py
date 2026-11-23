"""Compute and store Re G (eta-extrapolated) and branch spectral matrices for one frequency."""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bulk_green import total_G_eta, extrapolate, spectral_surface

omega = float(sys.argv[1])
N = int(sys.argv[2])
etas = [float(x) for x in sys.argv[3].split(",")]
nth = int(sys.argv[4]) if len(sys.argv) > 4 else 200
out = sys.argv[5] if len(sys.argv) > 5 else f"green_w{omega}_N{N}.npz"

t0 = time.time()
AT, AL, worst = spectral_surface(omega, ntheta=nth, nphi=2 * nth, check_monotone=True)
t1 = time.time()
G = total_G_eta(N, omega, etas)
t2 = time.time()
G0 = extrapolate(etas, G)
G0b = extrapolate(etas[1:], G[1:])          # drop the smallest eta: error indicator
imdiff = np.abs(G0.imag - (AT + AL)).max()
exdiff = np.abs(G0 - G0b).max()
print(f"omega={omega} N={N} etas={etas} ntheta={nth}")
print(f"  surface root residual {worst:.2e}; times: surface {t1-t0:.1f}s, grid {t2-t1:.1f}s")
print(f"  max|Im G_extrap - (A_T + A_L)| = {imdiff:.3e}   (scale Im G_xx(0) = {AT[5,1,1,0,0]+AL[5,1,1,0,0]:.6e})")
print(f"  max|G_extrap(all eta) - G_extrap(without smallest eta)| = {exdiff:.3e}")
print(f"  Re G_xx(0) = {G0[5,1,1,0,0].real:.8f}")
np.savez(out, omega=omega, N=N, etas=etas, ReG=G0.real, ImG=G0.imag, AT=AT, AL=AL)
