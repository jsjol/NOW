# NOW Migration Progress

## Current Phase: Phase 5 (nearly complete) — Utilities and Cleanup

### 2026-04-25

#### Phase 1: Foundation — COMPLETE
- Created project structure: `now/constraints/`, `tests/fixtures/`, `pyproject.toml`
- Installed JAX 0.10.0, pytest 9.0.3, pytest-cov in `now` conda env
- Created MATLAB test fixture generator (`generate_test_fixtures.m`) — **needs to be run in MATLAB**
- Consolidated architecture: `now/` is the package, `now_python/` to be removed after cleanup

#### Phase 2: Config and Timing — COMPLETE
- `now/config.py`: Full `NOW_config` class with all features:
  - `getActualTimings` with symmetry support
  - Motion compensation validation (infers linear/nonlinear from maxMagnitude)
  - Background compensation (general timing, specific timing)
  - Cross-term compensation (forces velocity comp when needed)
  - Signs vector for spin dephasing direction
- Tests: 16 tests in `test_config.py` (10 pass, 6 skipped awaiting MATLAB fixtures)

#### Phase 3: Constraints — COMPLETE
- `now/constraints/linear.py`: All linear constraints with all features:
  - Slew rate, max-norm gradient, echo condition, zero-gradient intervals
  - Linear motion compensation, background compensation, symmetry enforcement
  - `get_linear_constraint_matrices()` for testing against MATLAB
- `now/constraints/nonlinear.py`: All nonlinear constraints with analytical Jacobians:
  - Tensor encoding, gradient norm (L2), power/energy (per axis), Maxwell, motion compensation
  - `evaluate_all_nonlinear()` and `evaluate_all_nonlinear_jacobian()`
- Tests: 17 tests in `test_linear_constraints.py` + `test_nonlinear_constraints.py`
  - Jacobians verified against finite differences (atol=1e-3)
  - Jacobians verified against JAX autodiff (decimal=6)
  - 9 pass, 8 skipped awaiting MATLAB fixtures

#### Phase 4: Optimizer and Results — COMPLETE
- `now/optimize.py`: Full optimization pipeline:
  - SLSQP and trust-constr methods
  - Retry logic (up to 10 attempts with random restarts)
  - Initial guess scaling by target tensor diagonal
- `now/result.py`: `NOWResult` dataclass with all fields matching MATLAB output
  - Physical unit conversions (gamma, ms↔s, mT/m↔T/m)
  - md-dMRI compatible output (gwf, rf, dt)
- Tests: 5 end-to-end tests (all pass)
  - STE and LTE optimization complete successfully
  - STE produces approximately isotropic encoding tensor
  - All result fields populated

#### Phase 5: Utilities and Cleanup — IN PROGRESS
- `now/utils.py`: Utility functions ported (gamma, derivative matrices, integration weights, gwf_to_q, maxwell_coeff, write_wf, problem_to_name)
- `now/visualization.py`: Plot function ported from `now_python/`
- `now_example.py`: Standalone demo script
- Tests: 10 utility tests (all pass)
- **TODO**: Remove `now_python/` directory

#### Test Summary
- **48 tests pass** (config, constraints, objective, optimizer, utils)
- **16 tests skipped** (awaiting MATLAB fixtures from `generate_test_fixtures.m`)
- **0 failures**

### Next Steps
1. Run `generate_test_fixtures.m` in MATLAB to produce `.mat` reference files
2. Verify all 16 skipped tests pass against MATLAB fixtures
3. Remove `now_python/` directory
4. Merge `python` → `master` (manual intervention required)
5. Begin Phase 7: Refactoring on new branch
