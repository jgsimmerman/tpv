"""Self-contained arithmetic and reduced-model checks for TPV thesis draft.
Requires Python 3, NumPy and SciPy. No claim of full TPV validation.
"""
from fractions import Fraction as F
import json, math
import numpy as np
from scipy.integrate import quad

out={}
# Exact Haar-normalized strong-coupling coefficient.
C16=F(45,8)
blocks=[[(10,F(9,2)),(120,F(21,2)),(126,F(25,2))],
        [(1,F(0)),(45,F(8)),(210,F(12))]]
for block in blocks: assert sum(F(d,256) for d,c in block)==1
A=sum((F(d,256)/(2+c/C16) for block in blocks for d,c in block),F(0))
c2=4*(F(1,2)-A)
assert A==F(662262,1315237)
assert c2==F(-18574,1315237)
out['strong_coupling']={'A':str(A),'c2':str(c2),'decimal':float(c2)}
# Positive common kernel and the quartic no-fold polynomial.
beta=.05330497; z=.613; ell=4.2; C=1.45735
kernel=quad(lambda s:6*(s-1)*(2-s)*z/(z+beta*s),1,2,epsabs=1e-13)[0]
coef=[1,2-3*C,1+C]
disc=coef[1]**2-4*coef[0]*coef[2]
assert disc<0 and coef[2]>0
out['kernel']={'K':kernel,'correlation_lengths':[ell/math.sqrt(2*beta),ell/math.sqrt(beta)],'no_fold_polynomial':coef,'discriminant':disc}
# Exact reduced auxiliary Hamiltonian.
g=.1; J=1.
M=np.array([[0,3*np.sqrt(5)/2,0],[3*np.sqrt(5)/2,0,np.sqrt(42)],[0,np.sqrt(42),7.5]])
H=J*np.diag([-45/8,-13/8,3/8])+g*M
E,V=np.linalg.eigh(H)
prob=V[:,0]**2
assert abs(prob.sum()-1)<1e-13
out['auxiliary']={'energy':float(E[0]),'probabilities':prob.tolist()}
# Mixed sixth-order bounds: do not invent unknown source weights.
threshold=.097046248308400
out['mixed_shell']={'positive_if_norm2_below':(22/5)*threshold,'negative_if_norm2_above':14*threshold,'degree19_bound':(12/23)**20,'degree29_bound':(12/23)**30}
# Radial rotational benchmark from retrieved Q and I.
alpha=.005; LQ2=.059333; Q=.534166573628; I=94.01076
out['rotation']={'lambda_km2':alpha/LQ2,'delta_I_km3':alpha/LQ2*Q/4,'fraction':alpha/LQ2*Q/(4*I)}
# Causal two-vector Schur complement and spectral decomposition.
chi=.75;m2=2.
for k2 in [.1,1.,3.]:
 K=np.array([[k2+chi*m2,-chi*m2],[-chi*m2,chi*(k2+m2)]])
 eff=K[0,0]-K[0,1]**2/K[1,1]
 spectral=1/(1+chi)/k2+chi/(1+chi)/(k2+(1+chi)*m2)
 assert np.linalg.eigvalsh(K).min()>0
 assert abs(1/eff-spectral)<1e-12
out['memory']={'residues':[1/(1+chi),chi/(1+chi)],'mass2':(1+chi)*m2}
# Illustrative momentum/energy accounting, not an experimental input audit.
out['illustrative_radiation']={'force_N':237e-6,'minimum_massless_radiation_power_W':237e-6*299792458}
print(json.dumps(out,indent=2))
