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

## Phase 7b: Test suite cleanup — COMPLETE

### Changes made
- Removed dead code: `config.py` symmetry ValueError (line 36) — exhaustive search confirmed the `getActualTimings` symmetry path always produces equal durations by construction (`min(dur1, dur2)` used for both halves)
- Consolidated duplicated `_count_nonlinear` / `_count_nonlinear_constraints` — single source of truth now in `nonlinear.py`, imported by `problem.py`
- Added 5 tests for `optimize.py` error paths (all-attempts-fail, exception-then-retry, verbose output, verbose exception, redo-if-failed verbose)
- Added test for `linear.py` symmetry error (odd N without zero-gradient pause)
- Added test for `problem.py` nonlinear count with motion compensation
- Coverage: 93% total (100% on all modules except visualization and a 2-line ImportError fallback in config_io.py)

### Remaining coverage gaps (acceptable)
- `config_io.py:34-35` — `except ImportError` in `load_config` YAML path, only reachable when pyyaml is NOT installed. Not worth adding a mock test.
- `visualization.py` (0%) — matplotlib plotting code. Would require mocking the graphics backend. Intentionally excluded.

### Test suite notes
- **131 tests total** (121 fast, 10 slow), 0 warnings with `-W error`
- The `test_frozen` test in `test_problem.py` tests Python's dataclass decorator rather than application logic. Harmless but could be removed if test count matters.
- Jacobian finite-difference tests (`test_nonlinear_constraints.py`) test the same pattern across 4 configs. Could be parametrized to reduce code, but each exercises different constraint paths so all are valuable.

## Phase 8: MATLAB refactoring

### Goal
Align the MATLAB codebase with the Python refactoring so both share config formats and have a consistent user-facing API. The MATLAB side does NOT need the full three-layer architecture (MATLAB has one solver: `fmincon`), but it should be able to read/write the shared config format and produce compatible output.

### Shared config format (JSON)

Python's `save_config()` / `load_config()` serialize `NOW_config` to JSON with camelCase keys matching MATLAB property names. The MATLAB side should be able to read and write the same files.

**JSON schema** (as produced by `NOW_config.to_dict()`):
```json
{
  "targetTensor": [[1,0,0],[0,1,0],[0,0,1]],
  "N": 77,
  "useMaxNorm": false,
  "gMax": 80,
  "sMax": 100,
  "durationFirstPartRequested": 28,
  "durationSecondPartRequested": 22,
  "durationZeroGradientRequested": 8,
  "eta": 1,
  "enforceSymmetry": false,
  "redoIfFailed": true,
  "name": "NOW",
  "initialGuess": "random",
  "doMaxwellComp": true,
  "MaxwellIndex": 100,
  "MaxFunEval": 100000,
  "MaxIter": 5000,
  "motionCompensation": null,
  "doBackgroundCompensation": 0,
  "startTime": 0
}
```

**Key conventions:**
- `targetTensor` is a 3×3 nested array (row-major). MATLAB should use `cell2mat` or direct indexing.
- `motionCompensation` is `null` when empty (no orders). When present: `{"order": [1, 2], "maxMagnitude": [0, 1e-4]}`. The `linear` field is NOT serialized — it is inferred from `maxMagnitude <= 0` on construction (both MATLAB and Python do this).
- `x0` is only present when `initialGuess` is `"user-provided"`. Serialized as a flat array of length `3*N+1`.
- Boolean values are JSON `true`/`false` (not MATLAB 0/1). MATLAB's `jsondecode` handles this correctly.
- `MaxFunEval` and `MaxIter` are integers (no scientific notation).

**MATLAB implementation needed:**
```matlab
function problem = loadConfig(path)
    % Read JSON config and construct optimizationProblem
    text = fileread(path);
    d = jsondecode(text);
    problem = optimizationProblem(d);
end

function saveConfig(problem, path)
    % Save optimizationProblem to JSON
    d = struct();
    d.targetTensor = problem.targetTensor;
    d.N = problem.N;
    % ... (all public properties)
    if isempty(problem.motionCompensation.order)
        d.motionCompensation = [];  % becomes null in JSON
    else
        d.motionCompensation = struct('order', problem.motionCompensation.order, ...
            'maxMagnitude', problem.motionCompensation.maxMagnitude);
    end
    text = jsonencode(d, 'PrettyPrint', true);
    fid = fopen(path, 'w'); fprintf(fid, '%s', text); fclose(fid);
end
```

**Important:** `optimizationProblem.m` currently uses `eval` to set fields from a struct (line 68). The `loadConfig` function can pass the decoded JSON struct directly to `optimizationProblem()` since the field names match. However, MATLAB's `jsondecode` converts JSON arrays to MATLAB arrays and JSON `null` to empty `[]`, which is the correct behavior.

### Problem export format (MAT)

Python's `export_problem()` writes a `.mat` file with:
```
n_vars          : scalar (3*N + 1)
A_ineq          : (m, n_vars) double — inequality constraint matrix
b_ineq          : (m,) double — inequality RHS
A_eq            : (p, n_vars) double — equality constraint matrix
b_eq            : (p,) double — equality RHS
nonlinear_at_zero   : (k,) double — nonlinear constraints evaluated at x=0
nonlinear_jac_at_zero : (k, n_vars) double — Jacobian evaluated at x=0
n_nonlinear     : scalar — number of nonlinear constraints
```

This allows MATLAB to:
1. Verify constraint matrices match its own assembly (`defineLinearInequalityConstraints`, `defineLinearEqualityConstraints`)
2. Cross-validate nonlinear constraint values at known points
3. Use the matrices directly with external solvers

**MATLAB loading:**
```matlab
data = load('problem.mat');
% data.A_ineq, data.b_ineq, data.A_eq, data.b_eq are directly usable
```

### API correspondence (current and target)

| Concept | MATLAB (current) | Python (current) | MATLAB (target) |
|---------|------------------|-------------------|-----------------|
| Config | `optimizationProblem()` | `NOW_config()` | `optimizationProblem()` (unchanged) |
| Run | `NOW_RUN(problem)` | `now_optimize(config)` | `NOW_RUN(problem)` (unchanged) |
| Save config | — | `save_config(config, 'c.json')` | `saveConfig(problem, 'c.json')` |
| Load config | — | `load_config('c.json')` | `loadConfig('c.json')` |
| Export problem | — | `export_problem(problem, 'p.mat')` | Not needed (MATLAB can export its own matrices) |
| Import problem | — | `np.load('p.npz')` | `load('problem.mat')` |

### MATLAB changes needed (minimal)

1. **Add `saveConfig.m` and `loadConfig.m`** — JSON I/O functions (see templates above)
2. **Add `generate_test_fixtures.m` updates** — Export constraint matrices, nonlinear values, and Jacobians at known x-points for cross-validation. The existing `generate_test_fixtures.m` already does most of this.
3. **No changes to `optimizationProblem.m`** — The constructor already accepts a struct and sets fields. JSON-decoded structs are compatible.
4. **No changes to `optimize.m` / `NOW_RUN.m`** — The solver layer stays as `fmincon`. The MATLAB side doesn't need the solver abstraction.

### Cross-validation workflow

```
Python                          MATLAB
  │                               │
  ├── save_config(c, 'c.json') ──→ loadConfig('c.json')
  │                               ├── optimizationProblem(decoded)
  │                               ├── verify derived values match
  │                               │
  ├── export_problem(p, 'p.mat')─→ load('p.mat')
  │                               ├── compare A_ineq, b_ineq, A_eq, b_eq
  │                               ├── compare nonlinear values at x=0
  │                               │
  │←── load('matlab_result.mat')──├── save('matlab_result.mat', ...)
  │                               │
  └── compare b-values, tensors   └──
```

### Dead code and divergences to address in MATLAB

1. **`optimizationProblem.m` line 35-36:** The symmetry ValueError in `getActualTimings.m` (`'The first and second halves need to be equally long to enforce symmetry.'`) is dead code, same as the Python side — the algorithm takes `min(dur1, dur2)` ensuring equality by construction. Can be removed for consistency, or kept as a defensive guard.
2. **Hardcoded gamma:** MATLAB's `nonlconAnalytic.m` uses `2.6751e+08` (gamma_rad). Python consolidated this into `constants.py`. MATLAB should define `NOW_GAMMA_RAD = 42.576e6 * 2 * pi` in a shared location.
3. **`optimizationProblem.m` line 68:** Uses `eval()` to set fields from a struct — functional but fragile. Could be replaced with dynamic field access: `obj.(fieldNames{i}) = settings.(fieldNames{i})`.

## Phase 9: Future (do not implement now)

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

### Remaining improvements
- `optimize.py` could pass `ProblemParams` to `build_result()` instead of the full config to complete the decoupling — requires signature change in `build_result`
- Consider making `NOW_config` immutable (frozen dataclass or `__setattr__` override) to prevent accidental mutation after construction
