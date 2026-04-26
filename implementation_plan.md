# NOW Migration Implementation Plan

## Overview
Migrate NOW from MATLAB to Python using a parallel implementation with cross-validation strategy. Both codebases coexist and are validated against each other at every intermediate computation step.

## Phase 1: Foundation (Test Infrastructure & Directory Consolidation)
- [x] Set up conda environment with JAX, pytest
- [x] Create directory structure: `now/constraints/`, `tests/fixtures/`
- [x] Create `pyproject.toml`
- [x] Create MATLAB test fixture generator (`generate_test_fixtures.m`)
- [ ] **ACTION REQUIRED**: Run MATLAB fixture generator to produce `.mat` reference files
- [x] Consolidate `now/` and `now_python/` (keep `now/`, port best of `now_python/`)

## Phase 2: Config and Timing (TDD)
- [x] Write `test_config.py` — test `getActualTimings` against MATLAB fixtures
- [x] Validate `NOW_config` derived quantities against MATLAB
- [x] Implement `enforceSymmetry` in `getActualTimings`
- [x] Implement `motionCompensation` validation in `NOW_config`
- [x] Implement `doBackgroundCompensation` setup in `NOW_config`

## Phase 3: Constraints (TDD)
- [x] Write `test_linear_constraints.py` — compare assembled matrices vs MATLAB
- [x] Write `test_nonlinear_constraints.py` — values + Jacobians vs MATLAB and JAX
- [x] Implement `now/constraints/linear.py` with all features
- [x] Implement `now/constraints/nonlinear.py` with analytical Jacobians
- [x] Write JAX verification functions for constraint Jacobians
- [x] Implement missing: linear motion compensation, background compensation, symmetry, Maxwell

## Phase 4: Optimizer and Results (TDD)
- [x] Write `test_optimize.py` — end-to-end tests
- [x] Implement `now/optimize.py` with retry logic, initial guess, solver selection
- [x] Implement `now/result.py` — result dataclass with physical unit conversions
- [x] End-to-end validation: constraints satisfied, reasonable b-values

## Phase 5: Utilities and Cleanup
- [x] Port utility functions to `now/utils.py`
- [x] Port visualization to `now/visualization.py`
- [ ] Remove `now_python/` directory (pending: after MATLAB fixture validation)
- [x] Create standalone demo script (`now_example.py`)

## Phase 6: Cross-Validation and Merge
- [ ] Run `generate_test_fixtures.m` in MATLAB
- [ ] Verify all 16 skipped tests pass against MATLAB fixtures
- [ ] Remove `now_python/` directory
- [ ] Full test suite passes (48 + 16 = 64 tests)
- [ ] **PAUSE**: Merge `python` → `master` (manual intervention required)

## Phase 7: Refactoring (post-merge, new branch)
- [ ] Solver backend abstraction (Protocol-based)
- [ ] JAX backend implementation
- [ ] IPOPT backend implementation
- [ ] MATLAB-side refactoring
- [ ] Type hints, logging, config validation

## Phase 8: New Functionality (future, details TBD)
- JAX JIT-compiled optimization
- Automatic differentiation as default
- GPU support
- Manifold optimization exploiting encoding constraint structure
