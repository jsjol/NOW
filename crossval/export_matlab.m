% EXPORT_MATLAB  Export MATLAB NOW results for cross-validation against Python.
%
% Run from the NOW/ directory:
%   cd NOW
%   run crossval/export_matlab.m
%
% Requires the NOW MATLAB code on the path (optimizationProblem.m, etc.).

%% Add NOW root and private to path
nowRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(nowRoot);
addpath(fullfile(nowRoot, 'private'));

%% Create default problem
pObj = optimizationProblem();
N = pObj.N;

%% Deterministic x0 — linspace avoids RNG mismatch between MATLAB and Python
x0 = linspace(-1, 1, 3*N+1)';

% Apply the same target tensor scaling as getInitialGuess in optimize.m
% Note: MATLAB getInitialGuess uses (N+1):(2*N+1) and (2*N+1):(end-1),
% so element 2*N+1 gets scaled by both targetTensor(5) and targetTensor(9).
% With identity tensor (default) this is a no-op. We replicate exactly.
x0(1:N)               = x0(1:N)               * pObj.targetTensor(1);
x0((N+1):(2*N+1))     = x0((N+1):(2*N+1))     * pObj.targetTensor(5);
x0((2*N+1):(end-1))   = x0((2*N+1):(end-1))   * pObj.targetTensor(9);

%% Config-derived values
dt               = pObj.dt;
gMaxConstraint   = pObj.gMaxConstraint;
sMaxConstraint   = pObj.sMaxConstraint;
integralConstraint = pObj.integralConstraint;
tolMaxwell       = pObj.tolMaxwell;
signs            = pObj.signs;
zeroGradientAtIndex = pObj.zeroGradientAtIndex;  % 1-based
totalTimeActual  = pObj.totalTimeActual;

%% Derivative matrices (from optimize.m getDerivativeMatrices)
firstDerivativeMatrix = -diag(ones(N,1)) + diag(ones(N-1,1),1);
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);

secondDerivativeMatrix = diag(ones(N-1,1),-1) - 2*diag(ones(N,1)) + diag(ones(N-1,1),1);

%% Linear constraints (from optimize.m defineLinearInequalityConstraints/defineLinearEqualityConstraints)

% --- Inequality ---
if pObj.useMaxNorm
    A1 = kron(eye(3), firstDerivativeMatrix);
    A1 = [A1, zeros(size(A1,1),1)];
    b1 = gMaxConstraint * ones(size(A1,1),1);
else
    A1 = [];
    b1 = [];
end

A2 = kron(eye(3), secondDerivativeMatrix);
A2 = [A2, zeros(size(A2,1),1)];
b2 = sMaxConstraint * ones(size(A2,1),1);

A_ineq = [A1; -A1; A2; -A2];
b_ineq = [b1; b1; b2; b2];

% --- Equality ---
mc = pObj.motionCompensation;
if ~isempty(mc.linear)
    linear_ind = find(mc.linear);
else
    linear_ind = [];
end

n_eq_rows = 2 + length(zeroGradientAtIndex) + ...
    length(linear_ind) + ...
    1 * (pObj.doBackgroundCompensation > 0);

Aeq_single = zeros(n_eq_rows, N);

% Echo condition
Aeq_single(1, 1) = 1;
Aeq_single(2, N) = 1;

% Zero gradient at specified indices
if ~isempty(zeroGradientAtIndex)
    Aeq_single(2 + (1:length(zeroGradientAtIndex)), :) = ...
        firstDerivativeMatrix(zeroGradientAtIndex, :);
end

% Linear motion compensation
row_offset = 2 + length(zeroGradientAtIndex);
if ~isempty(linear_ind)
    t = ((1:N) - 1/2) * dt;
    for i = 1:length(linear_ind)
        order = mc.order(linear_ind(i));
        Aeq_single(row_offset + i, :) = -order * dt * t.^(order - 1);
    end
end

% Background compensation
if pObj.doBackgroundCompensation > 0
    row_bg = row_offset + length(linear_ind) + 1;
    s = pObj.startTime;
    H = cumsum([1; signs]) * dt + s;
    Aeq_single(row_bg, 1) = H(1) / 2;
    Aeq_single(row_bg, 2:(end-1)) = H(2:(end-1));
    Aeq_single(row_bg, end) = H(end) / 2;
end

% Symmetry enforcement
if pObj.enforceSymmetry
    if isempty(zeroGradientAtIndex)
        indicesBefore = 1:floor(N/2);
        indicesAfter = floor(N/2) + (1:ceil(N/2));
    else
        indicesBefore = 1:zeroGradientAtIndex(1);
        indicesAfter = (zeroGradientAtIndex(end)+1):N;
    end
    assert(length(indicesBefore) == length(indicesAfter), ...
        'Cannot enforce symmetry.');
    Nactivated = length(indicesBefore);
    Aeq_single = [Aeq_single; ...
        fliplr(eye(Nactivated)), zeros(Nactivated, N-2*Nactivated), -eye(Nactivated)];
end

A_eq = kron(eye(3), Aeq_single);
A_eq = sparse([A_eq, zeros(size(A_eq,1),1)]);
b_eq = zeros(size(A_eq,1),1);

% Convert sparse to full for saving
A_eq = full(A_eq);

%% Nonlinear constraints at x0
% nonlconAnalytic signature: (x, tolIsotropy, gMax, integralConstraint,
%   targetTensor, tolMaxwell*dt^2, signs, useMaxNorm, motionCompensation, dt)
tic;  % objFun requires a running timer
[c_raw, ceq, gradc_raw, gradceq] = nonlconAnalytic(x0, ...
    pObj.tolIsotropy, pObj.gMaxConstraint, pObj.integralConstraint, ...
    pObj.targetTensor, pObj.tolMaxwell * pObj.dt^2, ...
    pObj.signs, pObj.useMaxNorm, pObj.motionCompensation, pObj.dt);

c_nonlinear = c_raw;
% MATLAB gradc is (3N+1, n_constraints) — transpose to (n_constraints, 3N+1)
J_nonlinear = gradc_raw';

%% Objective at x0
[fval, grad] = objFun(x0);

%% Save
outDir = fullfile(fileparts(mfilename('fullpath')), 'data');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

save(fullfile(outDir, 'matlab_results.mat'), ...
    'x0', 'dt', 'gMaxConstraint', 'sMaxConstraint', ...
    'integralConstraint', 'tolMaxwell', 'signs', ...
    'zeroGradientAtIndex', 'totalTimeActual', ...
    'firstDerivativeMatrix', 'secondDerivativeMatrix', ...
    'A_ineq', 'b_ineq', 'A_eq', 'b_eq', ...
    'c_nonlinear', 'J_nonlinear', ...
    'fval', 'grad', 'N');

fprintf('Exported MATLAB results to %s\n', fullfile(outDir, 'matlab_results.mat'));
fprintf('  N = %d\n', N);
fprintf('  dt = %.10f\n', dt);
fprintf('  x0 size = [%d, %d]\n', size(x0));
fprintf('  c_nonlinear size = [%d, %d]\n', size(c_nonlinear));
fprintf('  J_nonlinear size = [%d, %d]\n', size(J_nonlinear));
fprintf('  A_ineq size = [%d, %d]\n', size(A_ineq));
fprintf('  A_eq size = [%d, %d]\n', size(A_eq));
fprintf('  zeroGradientAtIndex (1-based) = [');
fprintf(' %d', zeroGradientAtIndex);
fprintf(' ]\n');
