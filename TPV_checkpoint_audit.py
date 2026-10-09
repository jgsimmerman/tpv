"""Reproduce the three October 9 checkpoints from preserved originals.

Requires NumPy, SciPy, SymPy. Run from an extracted research package:
  python3 TPV_checkpoint_audit.py
Only /mnt/data output/input routing is changed in executed source strings.
The fixed-reference auxiliary data are supplied as a dependency, not rederived.
"""
from pathlib import Path
import contextlib, io, json, math, shutil, hashlib
import numpy as np
import scipy, sympy

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'supplements'
if not INPUT.exists():
    INPUT = ROOT.parent / 'upload'
DEPENDENCY = ROOT / 'dependencies' / 'mainTPV_coupled_throat_auxiliary_BVP_results_20261009.json'
if not DEPENDENCY.exists():
    DEPENDENCY = ROOT.parent / 'revision_inputs' / DEPENDENCY.name
RUN = ROOT / 'validation' / 'v03_reproduction'
RUN.mkdir(parents=True, exist_ok=True)
shutil.copyfile(DEPENDENCY, RUN / DEPENDENCY.name)
out = {'versions': {'numpy': np.__version__, 'scipy': scipy.__version__, 'sympy': sympy.__version__},
       'method': 'Preserved scripts executed with only /mnt/data rerouted; fixed-reference dependency not rederived.',
       'source_sha256': {}, 'checks': {}}
cases = [
    ('mainTPV_curvature_strain_principal_symbol_audit_20261009(1).py',
     'mainTPV_curvature_strain_principal_symbol_results_20261009.json'),
    ('mainTPV_tensor_spectral_parent_audit_20261009(1).py',
     'mainTPV_tensor_spectral_parent_results_20261009.json'),
    ('mainTPV_nonlinear_sGB_scan_20261009(1).py',
     'mainTPV_nonlinear_sGB_scan_results_20261009.json')]
for script, result in cases:
    original = (INPUT / script).read_text()
    out['source_sha256'][script] = hashlib.sha256((INPUT / script).read_bytes()).hexdigest()
    source = original.replace('/mnt/data', str(RUN))
    ns = {'__name__': '__main__', '__file__': str(INPUT / script)}
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        exec(compile(source, str(INPUT / script), 'exec'), ns)
    (RUN / (script + '.log')).write_text(log.getvalue())
    actual = json.loads((RUN / result).read_text())
    expected = json.loads((INPUT / result.replace('.json', '(1).json')).read_text())
    if 'nonlinear_sGB' not in script:
        errors = []
        def compare(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for key in a: compare(a[key], b[key])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for x, y in zip(a, b): compare(x, y)
            elif isinstance(a, (int, float)):
                assert math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10), (a, b)
                errors.append(abs(a-b))
            else:
                assert a == b
        compare(actual, expected)
        out['checks'][script] = {'assertions': 'passed', 'numeric_values': len(errors),
                                 'max_absolute_difference': max(errors)}
    else:
        rows, old = actual['rows'], expected['rows']
        assert len(rows) == len(old) == 30
        counts = np.zeros(5, dtype=int)
        endpoint_diffs = []
        for row, ref in zip(rows, old):
            assert 'error' not in row, row
            assert (row['eta'],row['d']) == (ref['eta'],ref['d'])
            assert row['events'] == ref['events']
            assert row['success'] and sum(row['events']) == 1
            counts += row['events']
            endpoint_diffs.append(abs(row['end_rho']-ref['end_rho']))
        assert counts.tolist() == [18,0,2,0,10]
        assert max(endpoint_diffs) < 1e-4
        # Reflection-even throat: all three first derivatives vanish and phi=1.
        throat = {ns['pr']:0, ns['pn']:0, ns['pu']:0, ns['uu']:1}
        d_symbol = next(q for q in ns['mat'].free_symbols if q.name == 'd')
        d_derivatives = [sympy.simplify(sympy.diff(eq.subs(throat),d_symbol)) for eq in ns['e']]
        assert all(q == 0 for q in d_derivatives)
        out['checks'][script] = {
            'runs':30, 'event_counts': counts.tolist(),
            'event_order':['low_lapse','small_radius','scalar_amplitude','high_lapse','radial_matrix'],
            'max_endpoint_rho_difference':max(endpoint_diffs),
            'max_local_throat_residual':max(r['throat_residual'] for r in rows),
            'max_sampled_constraint':max(r['max_constraint'] for r in rows),
            'symbolic_throat_equations_independent_of_d':True,
            'solver_success_is_terminal_event_not_global_solution':True}
        out['scan_by_eta'] = []
        for eta in [0., .05, .1]:
            group = [r for r in rows if r['eta']==eta]
            out['scan_by_eta'].append({'eta':eta,'c':group[0]['c'],
                'events':np.sum([r['events'] for r in group],axis=0).tolist(),
                'rho_range':[min(r['end_rho'] for r in group),max(r['end_rho'] for r in group)],
                'max_sampled_constraint':max(r['max_constraint'] for r in group)})

R=3.7436; xi=(1-.84739)/2; ke=10.; c=math.sqrt(1+2*xi*xi)
delta=-.117659*R**3/(2*ke*c)
out['cross_checkpoint']={'conditional_flare':xi,
    'previous_Kretschmann_lock_high_gradient_coefficient':ke*((2*xi/(R*c))**2+delta*2/c**3),
    'scope':'Same reduced ultrastatic homothetic calibration; not the nonlinear EsGB action.'}
out['scope']='Restricted script reproduction; no independent full field-equation, global BVP or Lorentzian stability validation.'
(ROOT / 'TPV_checkpoint_audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
