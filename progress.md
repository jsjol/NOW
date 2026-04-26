# NOW Migration & Refactoring Progress

## Phase 1–5: MATLAB-to-Python Migration — COMPLETE

### 2026-04-25
- Created project structure, conda env, pyproject.toml
- Implemented full `NOW_config` with all features (symmetry, motion comp, background comp)
- Implemented all linear constraints (slew, max-norm, echo, zero-gradient, motion, background, symmetry)
- Implemented all nonlinear constraints with analytical Jacobians (tensor, gradient norm, power, Maxwell, motion)
- Implemented optimization pipeline (SLSQP + trust-constr) with retry logic
- Implemented `NOWResult` with physical unit conversions and md-dMRI compatibility
- Ported all utility functions and visualization
- Created MATLAB test fixture generator (`generate_test_fixtures.m`)

## Phase 6: Cross-Validation and Merge — COMPLETE

### 2026-04-26
- [x] Generated MATLAB test fixtures
- [x] Verified MATLAB-agreement tests pass (16 fixture-dependent tests)
- [x] Added missing tests for utils.py (100%), linear.py (99%), nonlinear.py (99%), config.py (99%)
- [x] Fixed trust-constr bug (sparse/dense Jacobian mismatch) by refactoring `get_linear_constraints` to delegate to `get_linear_constraint_matrices`
- [x] Created bidirectional cross-validation protocol (`crossval/`)
- [x] Added `.coverage` and `crossval/data/` to `.gitignore`
- [x] Merge complete

**Test summary at merge:** 90 tests, 0 failures, 89% coverage (excluding visualization)

## Phase 7: Python Refactoring — IN PROGRESS

### 2026-04-26 (continued)

#### Step 1: Constants consolidation — COMPLETE
- [x] Created `now/constants.py` (GAMMA_HZ, GAMMA_RAD)
- [x] Updated `utils.py`, `result.py`, `nonlinear.py` to use constants
- [x] Added test_gamma_rad to test_utils.py
- [x] All 83 fast tests pass

#### Step 2: Problem layer — COMPLETE
- [x] Created `now/problem.py` with LinearConstraints, NonlinearConstraints, ProblemParams, OptimizationProblem, build_problem()
- [x] Created `tests/test_problem.py` (11 tests)
- [x] All 94 fast tests pass

#### Step 3: Solver layer — COMPLETE
- [x] Created `now/solvers/` package (protocol.py, scipy_solver.py, __init__.py)
- [x] Created `tests/test_solvers.py` (8 tests: 6 fast, 2 slow)
- [x] All 100 fast tests pass

#### Step 4: Rewire optimize.py — COMPLETE
- [x] Rewrote `now/optimize.py` to use build_problem() and get_solver()
- [x] Moved `objective()` from optimize.py to problem.py (fixed circular import)
- [x] Updated `__init__.py` imports
- [x] All 110 tests pass (100 fast + 10 slow)

#### Step 5: Clean up constraint modules — COMPLETE
- [x] Removed `get_linear_constraints(config, method)` from linear.py (solver-specific logic now in ScipySolver)
- [x] Removed `get_nonlinear_constraints(config, method)` from nonlinear.py
- [x] Removed `get_constraints()` wrapper from constraints/__init__.py
- [x] Removed scipy.optimize imports from constraint modules
- [x] Updated tests (removed 2 wrapper tests, renamed test class)
- [x] All 98 fast tests pass

#### Step 6: Config serialization — COMPLETE
- [x] Added `to_dict()` and `from_dict()` to NOW_config
- [x] Created `now/io/` package with `config_io.py` (save_config, load_config for JSON/YAML)
- [x] Created `tests/test_io.py` (12 tests: 8 dict round-trips, 4 file I/O)
- [x] All 109 fast tests pass (1 skipped — YAML not installed)

#### Step 7: Problem export — COMPLETE
- [x] Created `now/io/problem_io.py` with `export_problem()` (npz and mat formats)
- [x] Added 4 problem export tests to test_io.py
- [x] All 113 fast tests pass

#### Step 8: Update __init__.py exports — COMPLETE
- [x] Updated `now/__init__.py` with all new public exports
- [x] Example script runs successfully

#### Post-step cleanup — COMPLETE
- [x] Fixed Maxwell Jacobian divide-by-zero at x=0 (guard for m≈0 in nonlinear.py)
- [x] Suppressed scipy trust-constr `delta_grad == 0.0` warning in ScipySolver
- [x] Installed pyyaml in conda environment, YAML test now runs (was skipping)
- [x] Updated README.md with new features (solver backends, config serialization, problem export)
- [x] Updated MATLAB/Python correspondence table
- [x] Verified all module-level docstrings present
- [x] Final test run: 124 passed, 0 skipped, 0 warnings (with `-W error`), 89% coverage
