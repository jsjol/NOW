% generate_test_fixtures.m
% Generates reference data for Python cross-validation tests.
% Run from the NOW directory: >> generate_test_fixtures
%
% Output: tests/fixtures/*.mat files containing intermediate results
% for canonical problem configurations.

function generate_test_fixtures()

outputDir = fullfile(fileparts(mfilename('fullpath')), 'tests', 'fixtures');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Configuration 1: STE (Spherical Tensor Encoding) — default parameters
disp('Generating fixture: STE default...')
problem = optimizationProblem();
export_config_fixture(problem, fullfile(outputDir, 'config_ste_default.mat'));
export_constraint_fixture(problem, fullfile(outputDir, 'constraints_ste_default.mat'));

%% Configuration 2: STE with Maxwell compensation off
disp('Generating fixture: STE no Maxwell...')
settings = struct();
settings.doMaxwellComp = false;
problem_nomax = optimizationProblem(settings);
export_config_fixture(problem_nomax, fullfile(outputDir, 'config_ste_nomax.mat'));
export_constraint_fixture(problem_nomax, fullfile(outputDir, 'constraints_ste_nomax.mat'));

%% Configuration 3: LTE (Linear Tensor Encoding)
disp('Generating fixture: LTE...')
settings = struct();
settings.targetTensor = diag([1 0 0]);
problem_lte = optimizationProblem(settings);
export_config_fixture(problem_lte, fullfile(outputDir, 'config_lte.mat'));
export_constraint_fixture(problem_lte, fullfile(outputDir, 'constraints_lte.mat'));

%% Configuration 4: PTE (Planar Tensor Encoding)
disp('Generating fixture: PTE...')
settings = struct();
settings.targetTensor = diag([1 1 0]);
problem_pte = optimizationProblem(settings);
export_config_fixture(problem_pte, fullfile(outputDir, 'config_pte.mat'));
export_constraint_fixture(problem_pte, fullfile(outputDir, 'constraints_pte.mat'));

%% Configuration 5: STE with motion compensation
disp('Generating fixture: STE with motion compensation...')
settings = struct();
settings.motionCompensation = struct('order', [1 2], 'maxMagnitude', [0 1e-4], 'linear', []);
problem_motion = optimizationProblem(settings);
export_config_fixture(problem_motion, fullfile(outputDir, 'config_ste_motion.mat'));
export_constraint_fixture(problem_motion, fullfile(outputDir, 'constraints_ste_motion.mat'));

%% Configuration 6: STE with background compensation (general timing)
disp('Generating fixture: STE with background compensation...')
settings = struct();
settings.doBackgroundCompensation = 1;
settings.motionCompensation = struct('order', [1], 'maxMagnitude', [0], 'linear', []);
problem_bg = optimizationProblem(settings);
export_config_fixture(problem_bg, fullfile(outputDir, 'config_ste_bgcomp.mat'));
export_constraint_fixture(problem_bg, fullfile(outputDir, 'constraints_ste_bgcomp.mat'));

%% Configuration 7: STE with symmetry enforced
disp('Generating fixture: STE symmetric...')
settings = struct();
settings.enforceSymmetry = true;
settings.durationFirstPartRequested = 25;
settings.durationSecondPartRequested = 25;
settings.durationZeroGradientRequested = 8;
problem_sym = optimizationProblem(settings);
export_config_fixture(problem_sym, fullfile(outputDir, 'config_ste_symmetric.mat'));
export_constraint_fixture(problem_sym, fullfile(outputDir, 'constraints_ste_symmetric.mat'));

%% Configuration 8: STE with max-norm
disp('Generating fixture: STE max-norm...')
settings = struct();
settings.useMaxNorm = true;
problem_maxnorm = optimizationProblem(settings);
export_config_fixture(problem_maxnorm, fullfile(outputDir, 'config_ste_maxnorm.mat'));
export_constraint_fixture(problem_maxnorm, fullfile(outputDir, 'constraints_ste_maxnorm.mat'));

%% Configuration 9: Small N for quick tests
disp('Generating fixture: STE small N...')
settings = struct();
settings.N = 20;
problem_small = optimizationProblem(settings);
export_config_fixture(problem_small, fullfile(outputDir, 'config_ste_small.mat'));
export_constraint_fixture(problem_small, fullfile(outputDir, 'constraints_ste_small.mat'));

disp('All fixtures generated successfully!')
end


function export_config_fixture(problem, filepath)
% Export config-level data

config = struct();
config.N = problem.N;
config.dt = problem.dt;
config.gMax = problem.gMax;
config.sMax = problem.sMax;
config.eta = problem.eta;
config.targetTensor = problem.targetTensor;
config.useMaxNorm = problem.useMaxNorm;
config.enforceSymmetry = problem.enforceSymmetry;
config.doMaxwellComp = problem.doMaxwellComp;
config.MaxwellIndex = problem.MaxwellIndex;

config.gMaxConstraint = problem.gMaxConstraint;
config.sMaxConstraint = problem.sMaxConstraint;
config.integralConstraint = problem.integralConstraint;
config.tolMaxwell = problem.tolMaxwell;
config.tolIsotropy = problem.tolIsotropy;

config.durationFirstPartRequested = problem.durationFirstPartRequested;
config.durationSecondPartRequested = problem.durationSecondPartRequested;
config.durationZeroGradientRequested = problem.durationZeroGradientRequested;
config.durationFirstPartActual = problem.durationFirstPartActual;
config.durationSecondPartActual = problem.durationSecondPartActual;
config.durationZeroGradientActual = problem.durationZeroGradientActual;
config.totalTimeActual = problem.totalTimeActual;

config.zeroGradientAtIndex = problem.zeroGradientAtIndex;
config.signs = problem.signs;

config.motionCompensation = problem.motionCompensation;
config.doBackgroundCompensation = problem.doBackgroundCompensation;
config.startTime = problem.startTime;

save(filepath, '-struct', 'config');
end


function export_constraint_fixture(problem, filepath)
% Export constraint matrices and evaluated constraints at a deterministic x0

N = problem.N;

% Deterministic x0 for reproducibility
rng(42);
x0 = randn(3*N + 1, 1);
x0(1:N)           = x0(1:N)           * problem.targetTensor(1,1);
x0((N+1):(2*N))   = x0((N+1):(2*N))  * problem.targetTensor(2,2);
x0((2*N+1):(3*N)) = x0((2*N+1):(3*N)) * problem.targetTensor(3,3);

fixture = struct();
fixture.x0 = x0;
fixture.N = N;

% Derivative matrices
firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1);
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
secondDerivativeMatrix = diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1);

fixture.firstDerivativeMatrix = full(firstDerivativeMatrix);
fixture.secondDerivativeMatrix = full(secondDerivativeMatrix);

% Linear inequality constraints
if problem.useMaxNorm
    A1 = kron(eye(3), firstDerivativeMatrix);
    A1 = [A1 zeros(size(A1,1),1)];
    b1 = problem.gMaxConstraint*ones(size(A1,1),1);
else
    A1 = [];
    b1 = [];
end

A2 = kron(eye(3), secondDerivativeMatrix);
A2 = [A2 zeros(size(A2,1),1)];
b2 = problem.sMaxConstraint*ones(size(A2,1),1);

fixture.A_ineq = full([A1;-A1;A2;-A2]);
fixture.b_ineq = [b1;b1;b2;b2];

% Linear equality constraints — replicate the logic from optimize.m
Aeq_single = zeros(2 + length(problem.zeroGradientAtIndex) + ...
    nnz(problem.motionCompensation.linear) + ...
    1 * problem.doBackgroundCompensation, N);

Aeq_single(1,1) = 1;
Aeq_single(2,N) = 1;

if ~isempty(problem.zeroGradientAtIndex)
    Aeq_single(2+(1:length(problem.zeroGradientAtIndex)),:) = firstDerivativeMatrix(problem.zeroGradientAtIndex,:);
end

% Motion compensation (linear)
t = ((1:N)-1/2) * problem.dt;
linear_ind = find(problem.motionCompensation.linear);
for i = 1:length(linear_ind)
    order = problem.motionCompensation.order(linear_ind(i));
    Aeq_single(2+length(problem.zeroGradientAtIndex)+i,:) = - order * problem.dt * t.^(order-1);
end

% Background compensation
if problem.doBackgroundCompensation > 0
    s = problem.startTime;
    H = cumsum([1; problem.signs])' * problem.dt + s;
    row_idx = 2 + length(problem.zeroGradientAtIndex) + length(linear_ind) + 1;
    Aeq_single(row_idx, 1) = H(1)/2;
    Aeq_single(row_idx, 2:(end-1)) = H(2:(end-1));
    Aeq_single(row_idx, end) = H(end)/2;
end

% Symmetry
if problem.enforceSymmetry
    if isempty(problem.zeroGradientAtIndex)
        indicesBefore = 1:floor(N/2);
        indicesAfter = floor(N/2) + (1:ceil(N/2));
    else
        indicesBefore = 1:problem.zeroGradientAtIndex(1);
        indicesAfter = (problem.zeroGradientAtIndex(end)+1):N;
    end
    Nactivated = length(indicesBefore);
    Aeq_single = [Aeq_single; fliplr(eye(Nactivated)), zeros(Nactivated, N-2*Nactivated), -eye(Nactivated)];
end

Aeq = kron(eye(3), Aeq_single);
Aeq = [Aeq zeros(size(Aeq,1),1)];

fixture.A_eq = full(Aeq);
fixture.b_eq = zeros(size(Aeq,1),1);

% Nonlinear constraints at x0
[c, ceq, gradc, gradceq] = nonlconAnalytic(x0, problem.tolIsotropy, ...
    problem.gMaxConstraint, problem.integralConstraint, problem.targetTensor, ...
    problem.tolMaxwell*problem.dt^2, problem.signs, problem.useMaxNorm, ...
    problem.motionCompensation, problem.dt);

fixture.nonlinear_c = c;
fixture.nonlinear_ceq = ceq;
fixture.nonlinear_gradc = full(gradc);
fixture.nonlinear_gradceq = full(gradceq);

% Objective at x0
[f, g] = objFun(x0);
fixture.objective_f = f;
fixture.objective_g = g;

save(filepath, '-struct', 'fixture');
end
