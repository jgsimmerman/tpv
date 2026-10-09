"""TPV v0.2: source-comparison arithmetic and discrete spectrum normalization.
Run from the extracted research package; requires NumPy and SciPy.
These are restricted checks, not a full covariant-field-equation audit.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad

root = Path(__file__).resolve().parent
source = root / 'sources' / 'tpv_tensor_stability_audit.py'
if not source.exists():
    candidates = sorted((root / 'sources').glob('tpv_tensor_stability_audit*.py'))
    if not candidates:
        candidates = sorted((root.parent / 'sources').glob('tpv_tensor_stability_audit*.py'))
    source = candidates[0]
spec = importlib.util.spec_from_file_location('tt_archive', source)
tt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tt)
out = {}

# Use the same finite-difference GR operator, not continuum n*pi/R.
n = 800
base = tt.modes(n, 0)
k = np.arange(1, 6)
exact_discrete = 2 / base['dx'] * np.sin(k*np.pi/(2*(n+1)))
assert np.allclose(base['wave_numbers'], exact_discrete, rtol=1e-9)
out['tt'] = {'n_interior': n, 'GR_discrete_over_continuum': base['ratios'].tolist()}
for alpha in [.005, .05]:
    m = tt.modes(n, alpha)
    out['tt'][str(alpha)] = {
        'over_continuum': m['ratios'].tolist(),
        'over_discrete_GR': (m['wave_numbers']/base['wave_numbers']).tolist()}

c = 299792458.
G = 6.67430e-11
hbar = 1.054571817e-34
eps0 = 8.8541878128e-12
E = 25000/.00635
pressure = eps0*E**2/2
out['maxwell'] = {'E_V_per_m': E, 'pressure_Pa': pressure,
                  'area_for_237_microN_mm2': 237e-6/pressure*1e6}
out['casimir'] = []
for d in [4e-6, 1.5e-6, 1e-6, 1e-7]:
    u = np.pi**2*hbar*c/(720*d**4)
    out['casimir'].append({'d_m': d, 'abs_u_J_m3': u,
                          'h_scale': 8*np.pi*G*u*d*d/c**4})
L = 1.5e-6
u = np.pi**2*hbar*c/(720*L**4)
b = np.sqrt(u*32*np.pi*G*L**2*np.e/c**4)
j = c**4*b/(8*np.pi*G*L**2)
out['gaussian_shift'] = {'b_density_match': b, 'center_flux_energy_units_J_m3': j,
                         'flux_over_Casimir_density': j/u}
beta = .05330497
low_z = quad(lambda s: 6*(s-1)*(2-s)/(beta*s),1,2)[0]
out['kernel'] = {'linear_coefficient': low_z,
                 'at_1_micron_linear_approx': low_z*(4.2*1.616255e-35/1e-6)**2}

# Fixed densitized-current variation: include determinant and normalization.
g = np.diag([-1.,1.,1.,1.])
J = np.array([1.2,.2,.1,0.])
dg = np.array([[.1,.03,0,0],[.03,-.04,.02,0],[0,.02,.05,0],[0,0,0,-.02]])
def fluid(metric):
    density = np.sqrt(-J@metric@J)/np.sqrt(-np.linalg.det(metric))
    velocity = J/(np.sqrt(-np.linalg.det(metric))*density)
    return density, velocity
density, velocity = fluid(g)
expected_dn = -density/2*np.sum((np.linalg.inv(g)+np.outer(velocity,velocity))*dg)
expected_du = velocity/2*(velocity@dg@velocity)
step = 1e-5
plus, minus = fluid(g+step*dg), fluid(g-step*dg)
err_n = abs((plus[0]-minus[0])/(2*step)-expected_dn)
err_u = np.max(abs((plus[1]-minus[1])/(2*step)-expected_du))
assert err_n < 1e-9 and err_u < 1e-9
out['fixed_current'] = {'density_derivative_abs_error': float(err_n),
                         'velocity_derivative_max_abs_error': float(err_u)}

# Synthetic compact-support projection: j=w0=1, phi=r^2(1-r)^4.
# H[phi]=(r^4 phi')'; boundary flux and weighted integral vanish.
phi = np.polynomial.Polynomial([0,0,1,-4,6,-4,1])
r4 = np.polynomial.Polynomial([0,0,0,0,1])
Hphi = (r4*phi.deriv()).deriv()
projection = quad(Hphi,0,1,epsabs=1e-12)[0]
assert abs(projection)<1e-12 and abs(Hphi(.5))>1e-4
out['synthetic_projection'] = {'integral': projection,
                                'local_source_at_half': Hphi(.5),
                                'scope': 'boundary cancellation illustration; not a stellar rerun'}
print(json.dumps(out,indent=2))
