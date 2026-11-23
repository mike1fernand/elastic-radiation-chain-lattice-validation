import json
from itertools import product
from pathlib import Path
import numpy as np

if not __debug__:
    raise SystemExit("Run with assertions enabled")

N = 6
L = N + 1
R = Path(__file__).resolve().parent / "upstream"

def contact_matrix(g):
    H = np.zeros((L + 2*N, L + 2*N))
    for s, strength in enumerate(g):
        for y, z in product((-1, 1), repeat=2):
            corner = np.array([s+.5, .5*y, .5*z])
            for n in (s, s+1):
                direction = corner - np.array([n, 0, 0])
                direction /= np.linalg.norm(direction)
                a = np.zeros(len(H))
                a[n] = -direction[0]
                a[L + 2*s] = direction[0]/2
                a[L + 2*s + 1] = np.dot(direction, [0,y,z])/(2*np.sqrt(2))
                H += strength * np.outer(a, a)
    return H

def solve(g, mats, incident_end):
    H = contact_matrix(g)
    G, *projected = (matrix[:2*N, :2*N] for matrix in mats)
    V = H[L:, :L]
    force_per_line = -np.linalg.solve(np.eye(2*N) + H[L:, L:] @ G, V)
    w, K = 1.2, 1.69
    z = 1 - w*w/(2*K)
    phase = z + 1j*np.sqrt(1-z*z)
    gamma = 2*K*phase.imag
    lead_operator = (w*w-2*K)*np.eye(L, dtype=complex)
    for j in range(L-1):
        lead_operator[j,j+1] = lead_operator[j+1,j] = K
    lead_operator[0,0] += K*phase
    lead_operator[-1,-1] += K*phase
    b = np.zeros(L)
    b[incident_end] = np.sqrt(gamma)
    q = np.linalg.solve(lead_operator-H[:L,:L]-V.T@G@force_per_line,b)
    f = force_per_line @ q
    far_end = L-1-incident_end
    ref = abs(1 - 1j*np.sqrt(gamma)*q[incident_end])**2
    tran = gamma*abs(q[far_end])**2
    powers = [2*np.vdot(f, Im@f).real for Im in projected]
    assert abs(sum(powers)+ref+tran-1) < 3e-13
    return np.asarray([ref,tran,sum(powers[:2]),powers[2]])


ramp = np.arange(1.0, 7.0) / 7
cases = {
    'ramp': (ramp,0),
    'same_entrance_swap_3_and_5': (ramp[[0,1,4,3,2,5]],0),
    'reversed': (ramp,L-1),
}
reference = json.loads((R / 'profile_reversal.json').read_text())
with np.load(R / 'REFERENCE_KERNEL.npz',allow_pickle=False) as matrix:
    surface = [matrix['G0'],*(matrix[f'Im_branch_{i}'] for i in range(3))]
with np.load(R / 'N192.npz',allow_pickle=False) as matrix:
    tetra = [matrix['G0'],*(matrix[f'Im{i}'] for i in range(3))]
for kind,mats in [('surface',surface),('tetrahedron',tetra)]:
    values={name:solve(profile,mats,end) for name,(profile,end) in cases.items()}
    np.testing.assert_allclose(values['ramp'][1],values['reversed'][1],rtol=0,atol=1e-14)
    if kind=='surface':
        for name,v in values.items():
            np.testing.assert_allclose(v,[reference[name][k] for k in ('reflection','through','qT','qL')],rtol=0,atol=2e-13)
    print(kind,'qT estimates (percent of incident power)',{name:round(float(100*v[2]),7) for name,v in values.items()})
