"""Exact-rational residual certificate for the P1 swap of the third and fifth motif at the nominal point.

NumPy only proposes an inverse. Every accepted inequality is verified using
Fraction arithmetic and rational intervals for sqrt(2) and sqrt(133).
The cached bare Green and qT error caps are inputs inherited from the pinned
repository. The new calculation propagates them through the swapped contact.
No interval package or floating-point error assumption is used for acceptance.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import argparse, hashlib, json
import numpy as np

if not __debug__:
    raise SystemExit("Verification uses assertions. Run Python without -O or -OO.")
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("kernel",type=Path,help="Path to the pinned REFERENCE_KERNEL.npz")
SOURCE=parser.parse_args().kernel
HASH='82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==HASH

class C:
    __slots__=('r','i')
    def __init__(self,r=0,i=0):self.r=F(r);self.i=F(i)
    def __add__(self,b):
        if not isinstance(b,C):b=C(b)
        return C(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,b):return self+-asC(b)
    def __mul__(self,b):
        b=asC(b);return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):
        assert not isinstance(b,C);return C(self.r/F(b),self.i/F(b))
    def conj(self):return C(self.r,-self.i)
    def sq(self):return self.r*self.r+self.i*self.i
    def l1(self):return abs(self.r)+abs(self.i)
    def fl(self):return complex(float(self.r),float(self.i))

def asC(z):return z if isinstance(z,C) else C(z)
def zeros(m,n):return [[C() for _ in range(n)] for _ in range(m)]
def eye(n):return [[C(int(i==j)) for j in range(n)] for i in range(n)]
def minus(A,B):return [[a-b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def mm(A,B):
    bt=list(zip(*B));return [[sum((a*b for a,b in zip(row,col)),C()) for col in bt] for row in A]
def fro2(A):return sum(x.sq() for row in A for x in row)
def vec2(x):return sum(z.sq() for z in x)
def sqrt_bounds(n,digits=80):
    den=10**digits;lo=F(isqrt(n*den*den),den);return lo,lo+F(1,den)
S2L,S2U=sqrt_bounds(2);S13L,S13U=sqrt_bounds(133)
S2=(S2L+S2U)/2;S13=(S13L+S13U)/2
e2=(S2U-S2L)/2;e13=(S13U-S13L)/2
k=F(169,100);om=F(6,5);g=[F(j,7) for j in (1,2,5,4,3,6)]
with np.load(SOURCE,allow_pickle=False) as f:
    gf=f['G0'][:12,:12].copy();qf=f['Im_qT'][:12,:12].copy()
    assert np.array_equal(qf,f['Im_branch_0'][:12,:12]+f['Im_branch_1'][:12,:12])
G=[[C(F(float(z.real)),F(float(z.imag))) for z in row] for row in gf]
Q=[[F(float(z)) for z in row] for row in qf]
assert all(Q[i][j]==Q[j][i] for i in range(12) for j in range(12))

# Exact LDL proves that the dyadic output metric is positive definite.
L=[[F(int(i==j)) for j in range(12)] for i in range(12)];piv=[]
for j in range(12):
    d=Q[j][j]-sum(L[j][h]**2*piv[h] for h in range(j))
    assert d>0;piv.append(d)
    for i in range(j+1,12):
        L[i][j]=(Q[i][j]-sum(L[i][h]*L[j][h]*piv[h] for h in range(j)))/d
assert all(sum(L[i][h]*piv[h]*L[j][h] for h in range(12))==Q[i][j] for i in range(12) for j in range(12))

def model(s2,s13):
    U=zeros(12,12);V=zeros(12,7);Kself=[F(0)]*7
    for j,gj in enumerate(g):
        U[2*j][2*j]=C(2*gj/3);U[2*j+1][2*j+1]=C(4*gj/3)
        V[2*j][j]=V[2*j][j+1]=C(-2*gj/3)
        V[2*j+1][j]=C(-2*s2*gj/3);V[2*j+1][j+1]=C(2*s2*gj/3)
        Kself[j]+=4*gj/3;Kself[j+1]+=4*gj/3
    A=zeros(19,19);D=zeros(19,12)
    for i in range(7):
        for j in range(7):
            h=C(om*om-2*k-Kself[i]) if i==j else C(k if abs(i-j)==1 else 0)
            if i==j and i in [0,6]:h+=C(F(97,100),3*s13/25)
            A[i][j]=h/k
    for i in range(7):
        for j in range(12):
            A[i][7+j]=sum(V[h][i].conj()*G[h][j] for h in range(12))
            D[i][j]=V[j][i].conj()/k
    for i in range(12):
        for j in range(7):A[7+i][j]=-V[i][j]/k
        for j in range(12):
            A[7+i][7+j]=C(int(i==j))+sum(U[i][h]*G[h][j] for h in range(12))
            D[7+i][j]=U[i][j]/k
    return A,D

Ac,Dc=model(S2,S13);A00,D00=model(F(0),F(0))
A2,D2=model(F(1),F(0));A2=minus(A2,A00);D2=minus(D2,D00)
A13,_=model(F(0),F(1));A13=minus(A13,A00)
deltaA=e2*sum(z.l1() for row in A2 for z in row)+e13*sum(z.l1() for row in A13 for z in row)
deltaD=e2*sum(z.l1() for row in D2 for z in row)
assert deltaA<F(1,10**60) and deltaD<F(1,10**60)

Rf=np.linalg.inv(np.array([[z.fl() for z in row] for row in Ac],complex))
R=[[C(F(float(z.real)),F(float(z.imag))) for z in row] for row in Rf]
E=minus(eye(19),mm(R,Ac));eps_c=F(1,10**13);eps=F(1,10**12)
assert fro2(E)<eps_c**2

# eps bounds ||I - R A_c|| including the uncertainty in the exact radicals.
assert fro2(R)<36 and eps_c+6*deltaA<eps<1
mu=eps/(1-eps)
RD=mm(R,Dc); Hhat=RD[7:]
col=[R[i][0] for i in range(19)];v=col[7:]
assert fro2(RD)<4 and vec2(col)<4
def quadratic(v):
    out=sum((v[i].conj()*Q[i][j]*v[j] for i in range(12) for j in range(12)),C())
    assert out.i==0
    return out.r
# ||A_c^-1 D - R D_center|| <= mu*||R D_center|| + ||A_c^-1||*deltaD.
dc=mu*2 + 6*deltaD/(1-eps)
assert fro2(Hhat)<F(949,1000)**2 and F(949,1000)+dc<F(95,100)
# Weighted norm uses metric F_c^* F_c = 2k P_eta^* Q P_eta.
Fcap=F(2253,1000)
assert 2*k*sum(Q[i][i] for i in range(12)) < Fcap**2
jhat2=2*k*sum(quadratic(list(col)) for col in zip(*Hhat))
assert jhat2<F(6415,10000)**2
assert F(6415,10000)+Fcap*dc<F(642,1000)
gamlo=6*S13L/25; gamhi=6*S13U/25
c_cap=F(32,25);assert gamhi/k<c_cap**2
# Nominal source norm and cached amplitude interval, including inverse residual.
assert gamhi/k*vec2(v)<F(569,1000)**2
assert F(569,1000)+c_cap*2*mu < F(57,100)
amp2lo=2*gamlo*quadratic(v); amp2hi=2*gamhi*quadratic(v)
assert amp2lo>F(823754,10**6)**2 and amp2hi<F(823755,10**6)**2
ae=Fcap*c_cap*mu*2
alo=F(82375,100000); ahi=F(82376,100000)
assert F(823754,10**6)-ae>alo and F(823755,10**6)+ae<ahi
# Existing independent bare-bulk premises, scaled by k.
eg=F(1,250); eq=F(1,500); h=F(95,100); jcap=F(642,1000); z0=F(57,100)
assert h*eg<1
Z=z0/(1-h*eg)
d=jcap*eg*Z
assert alo>d
lo=(alo-d)**2-2*eq*Z**2
hi=(ahi+d)**2+2*eq*Z**2
# Strict comparisons, independent of rounding of the displays.
assert F(674,1000)<lo<hi<F(683,1000)
assert F(732,1000)-F(683,1000)==F(49,1000)
result={'claim':'pointwise bound for qT fraction at omega=6/5, same entrance, swap 3 and 5',
        'g':[str(x) for x in g], 'eG':str(eg),'eQ':str(eq),'kernel_sha256':HASH,
        'qT_lower':str(lo),'qT_upper':str(hi),'qT_lower_display':float(lo),'qT_upper_display':float(hi)}
print(json.dumps(result,indent=2))
