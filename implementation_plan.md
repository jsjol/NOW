# NOW Implementation Plan

## Phase 1–6: MATLAB-to-Python Migration — COMPLETE

All checkboxes completed. 90 tests pass, 89% coverage. See `progress.md` for details.

## Phase 7: Python Refactoring

### Goal
Refactor from MATLAB port to idiomatic, extensible Python. Separate config → problem → solver → result. Enable swappable solver backends and shared config formats.

### Architecture: Three-layer separation

```
User layer:     NOW_config  →  now_optimize()  →  NOWResult
                    ↓                                  ↑
Problem layer:  build_problem(config)          build_result(solver_result, problem)
                    ↓                                  ↑
Solver layer:   solver.solve(problem, x0)  →  SolverResult
```

### Steps

- [ ] **Step 1: Constants consolidation** — `now/constants.py` with `GAMMA_HZ`, `GAMMA_RAD`. Remove duplication from `utils.py`, `nonlinear.py`, `result.py`.
- [ ] **Step 2: Problem layer** — `now/problem.py` with `ProblemParams`, `LinearConstraints`, `NonlinearConstraints`, `OptimizationProblem`, `build_problem()`.
- [ ] **Step 3: Solver layer** — `now/solvers/` with `SolverProtocol`, `SolverResult`, `ScipySolver`, `get_solver()`.
- [ ] **Step 4: Rewire optimize.py** — Use `build_problem()` + `get_solver()` internally. Public API unchanged.
- [ ] **Step 5: Clean up constraints** — Remove `method` parameter from constraint functions. Solver-specific translation moves to `ScipySolver`.
- [ ] **Step 6: Config serialization** — `to_dict()`/`from_dict()` on `NOW_config`. JSON/YAML I/O in `now/io/config_io.py`.
- [ ] **Step 7: Problem export** — `now/io/problem_io.py` for exporting NLP to npz/mat.
- [ ] **Step 8: Type hints and exports** — Update `__init__.py`, add type hints throughout.

### Design decisions

| Decision | Rationale |
|----------|-----------|
| `OptimizationProblem` is a dataclass, not a class hierarchy | Different formulations → different instances via different builders |
| Keep derived values on `NOW_config` | Backward compat. `ProblemParams` provides clean frozen alternative |
| Constraints stay as monolithic functions | Splitting into per-constraint objects is premature without DNLP |
| Solver options as plain dicts | Simple, doesn't over-constrain backends |

### DNLP considerations

NOW's problem is entirely smooth (C²). DNLP's nonsmooth handling doesn't apply directly. But the architecture prepares for it: `OptimizationProblem` maps to NLP standard form, linear constraints map to DNLP atoms, and a future `DNLPSolver` could re-express the problem using CVXPY/DNLP. The nonconvex encoding constraint (X^T X = bB̂) stays as a general smooth constraint in any DNLP formulation.

## Phase 8: Future (do not implement now)

### Solver backends
- **DNLP backend**: Re-express problem via CVXPY, solve with Ipopt/Knitro
- **Ipopt backend**: Direct via `cyipopt`, two-sided inequalities with sparsity
- **JAX backend**: JIT-compiled objective + constraints, autodiff Jacobians

### Formulations and algorithms
- **Alternative formulations**: NOW-revisited (equality encoding, equality Maxwell, simplified heating) via `build_problem_v2(config)`
- **Manifold optimization**: Exploit Stiefel manifold structure of encoding constraint
- **Batch optimization**: Multiple waveforms simultaneously

### Applications
- **Closed-loop MRI**: Constraint-aware experiment design (PROFOUND connection)

### Improvements noticed during Phase 7
- `nonlinear.py` lines 180–189 (`_count_nonlinear_constraints`) duplicates `problem.py`'s `_count_nonlinear` — consolidate into one location
- Maxwell Jacobian `dc6_dM = (1/m) * M` produces divide-by-zero at x=0 — add guard for m≈0 (not a real optimization issue since m=0 is never a feasible point, but produces RuntimeWarning in tests)
- `config_io.py` YAML branch has 71% coverage (skipped when pyyaml not installed) — add pyyaml to dev dependencies or accept the skip
- `optimize.py` could pass `ProblemParams` to `build_result()` instead of the full config to complete the decoupling — requires signature change in `build_result`
- Consider making `NOW_config` immutable (frozen dataclass or `__setattr__` override) to prevent accidental mutation after construction
