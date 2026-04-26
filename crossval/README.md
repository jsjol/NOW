# NOW Cross-Validation: MATLAB vs Python

Bidirectional regression test: both codebases export results for identical
fixed inputs, and a comparison script checks they match.

## Usage

1. In MATLAB (from `NOW/`): `run crossval/export_matlab.m`
2. In a terminal (from `NOW/`): `python crossval/export_python.py`
3. Compare: `python crossval/compare_results.py`

Exit code 0 means all checks passed; exit code 1 means at least one failed.

Note: `zeroGradientAtIndex` is 1-based in MATLAB and 0-based in Python.
The comparison script accounts for this offset automatically.
