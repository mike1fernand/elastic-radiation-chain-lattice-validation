import hashlib
from pathlib import Path
from mpmath import iv
import numpy as np
if not __debug__:
 raise SystemExit("Run with assertions enabled")
iv.dps=60
def rat(n,d=1):return iv.mpf(n)/d
def dy(n):
 n,d=float(n).as_integer_ratio()
 return rat(n,d)
def complex_iv(t):return iv.mpc(dy(t.real),dy(t.imag))
def abs2(c):return c.real**2+c.imag**2
kernel=Path(__file__).resolve().parent/'upstream'/'REFERENCE_KERNEL.npz'
assert hashlib.sha256(kernel.read_bytes()).hexdigest()=='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
with np.load(kernel) as d:
 G=iv.matrix([[complex_iv(v) for v in row] for row in d['G0'][:12,:12]])
 Q=iv.matrix([[dy(v) for v in row] for row in d['Im_qT'][:12,:12]])
Gcap=rat(1,250);Qcap=rat(1,500)
k=rat(169,100);omega=rat(6,5);signal=iv.sqrt(6*iv.sqrt(133)/(25*k))
Fcap=iv.sqrt(2*k*sum(Q[j,j] for j in range(12))) # PSD previously certified; tr bounds largest eigenvalue
for ns in [[1,2,3,4,5,6],[1,2,5,4,3,6]]:
 U=iv.matrix(12,12);V=iv.matrix(12,7);h=iv.matrix(7,7)
 for j in range(7):
  h[j,j]=omega**2-2*k
  if j<6:h[j,j+1]=h[j+1,j]=k
 for j in [0,6]:h[j,j]+=iv.mpc(rat(97,100),3*iv.sqrt(133)/25)
 for j,n in enumerate(ns):
  c=rat(2*n,21);U[2*j,2*j]=c;U[2*j+1,2*j+1]=2*c;V[2*j,j]=V[2*j,j+1]=-c
  V[2*j+1,j]=-iv.sqrt(2)*c;V[2*j+1,j+1]=iv.sqrt(2)*c
  h[j,j]-=2*c;h[j+1,j+1]-=2*c
 # Use scaled matrices and select source dof
 A=iv.matrix(19,19);A[:7,:7]=h/k;A[:7,7:]=V.T*G;A[7:,:7]=-V/k;A[7:,7:]=iv.eye(12)+U*G
 D=iv.matrix(19,12);D[0:7,:]=V.T/k;D[7:,:]=U/k
 R=A**-1
 w=(R*D)[7:,:]
 Hcap=iv.sqrt(sum(abs2(w[i,j]) for i in range(12) for j in range(12)))
 z=R[7:,0]*signal
 z0=iv.sqrt(sum(abs2(x) for x in z))
 qc=Q*z
 # This simplified form is real because Q is real symmetric.
 tc=2*k*sum(z[j].real*qc[j].real+z[j].imag*qc[j].imag for j in range(12))
 Z=z0/(1-Hcap*Gcap)
 W=Q*w
 j2=2*k*sum(w[i,j].real*W[i,j].real+w[i,j].imag*W[i,j].imag for i in range(12) for j in range(12))
 J=iv.sqrt(j2)
 amp=iv.sqrt(tc)
 if ns[2]==5:
  assert Hcap<rat(95,100) and J<rat(642,1000) and z0<rat(57,100) and amp>rat(82375,100000) and amp<rat(82376,100000)
  print("tight caps verified")
 delta=J*Gcap*Z
 low=(iv.sqrt(tc)-delta)**2-2*Qcap*Z**2
 upper=(iv.sqrt(tc)+delta)**2+2*Qcap*Z**2
 print(ns, "Tcenter",tc,"H",Hcap,"Z",Z,"F",Fcap,"lo",low,"hi",upper)
 assert Hcap*Gcap<1 and iv.sqrt(tc)>delta
 if ns[2]==5:assert upper<rat(683,1000)
 else:assert low>rat(732,1000)
