# Tensor Polarizable Vacuum — consolidated research draft v0.3

Source cutoff: 9 October 2026. Author attribution: Jacob G. Simmerman.

## Deliverables
- TPV_Thesis.tex: self-contained thesis-style LaTeX manuscript.
- TPV_Thesis.pdf: compiled manuscript.
- TPV_thesis_audit.py: independent compact audit (Python 3, NumPy, SciPy).
- TPV_thesis_audit_results.json: captured output of that audit.
- TPV_source_inventory.csv: all 692 retrieved selected files, hashes, source identities, and duplicate mappings.
- sources/: one original file per distinct SHA-256 content (616 files).
- TPV_revision_audit.py and TPV_revision_audit_results.json: corrected discrete spectrum, laboratory scale calculations, and restricted variation/boundary tests.
- TPV_supplemental_source_manifest.json: hashes and identities of all eleven supplemental attachments and the recovered auxiliary-data dependency.
- CHANGELOG.md: edition changes and validation boundaries.
- supplements/: the monograph, duplicate D5 source, and nine checkpoint files preserved verbatim.
- dependencies/: the recovered auxiliary BVP results used by the nonlinear scan.
- TPV_checkpoint_audit.py and TPV_checkpoint_audit_results.json: three-checkpoint reproduction, symbolic throat test, and combined interpretation.
- requirements_checkpoints.txt: versions used for that reproduction.
- validation/: captured original audit runs, including the dependency failure described below.

## Scope
This is a critical consolidated draft, not a completed doctoral dissertation or proof of a unified TPV theory. The main text reconciles key corrections and explicitly identifies which claims are exact, conditional, historical, or open. Complete exports of every later TPV, cTPV, aTPV and electrogravity chat were not established. Preservation in sources/ does not imply line-by-line adjudication or execution of all code. Some earlier programs refer to missing absolute paths or modules. Original research files are preserved without silently rewriting their claims.

## Reproduction
Compile twice with pdfLaTeX:

    pdflatex -interaction=nonstopmode -halt-on-error TPV_Thesis.tex
    pdflatex -interaction=nonstopmode -halt-on-error TPV_Thesis.tex
    python3 TPV_thesis_audit.py
    python3 TPV_revision_audit.py
    python3 -m pip install -r requirements_checkpoints.txt
    python3 TPV_checkpoint_audit.py

Original code successfully rerun during earlier editions and retained here:
- tpv_tensor_stability_audit.py --n 800
- mainTPV_local_constraint_topology_audit_20261009.py
- mainTPV_Berry_exact_modular_certificate_20261009.py
- mainTPV_full_affine_hessian_audit_20261009.py

The archived mainTPV_YM_plaquette_normalization_audit.py could not run during the original build because SymPy was absent. SymPy was subsequently installed for the v0.3 checkpoint audit; the historical normalization script was not rerun. The compact audit independently reproduces its rational coefficient using fractions; no full archived-script pass is claimed.

No full geon PDE, unrestricted covariant stability proof, or active 3D merger was run. Bibliographic checks cover selected external methodology, not exhaustive novelty or experimental literature review.

The revised PDF incorporates the monograph without treating all of its reported checks as newly rerun. Its mislabeled TT table has been corrected by recomputing the corresponding discrete GR reference. The D5 attachment adds no new content. The original inventory and its 616 distinct source contents remain preserved; supplemental attachments are listed separately.

## Version 0.3 checkpoint reproduction

The curvature-locking audit (89 numeric values) and tensor spectral-parent audit (38 numeric values) reproduce their supplied JSON exactly. All 30 nonlinear EsGB scan events and terminal radii also reproduce exactly. Eighteen runs reach the low-lapse cutoff, two the scalar-amplitude cutoff, and ten the static radial-matrix threshold. None establishes the intended global exterior. The wrapper separately verifies that all four reflection-symmetric throat equations are independent of the nonlinear coupling parameter d.

The wrapper redirects only historical /mnt/data paths, preserving originals under supplements/. Its run products are in validation/v03_reproduction/. The fixed-reference auxiliary BVP data are inputs, not a newly solved global BVP. The separate collocation failures mentioned in the checkpoint remain reported because that implementation was not supplied.

The restricted static locking Hessian, physical-polarization spectral diagnostic, and nonlinear scalar control are different models. No full Lorentzian TPV stability, one-ended geon, Standard Model emergence, or active merger success is claimed.
