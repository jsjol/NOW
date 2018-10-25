function [result, problem] = optimize(varargin)
% Script that optimizes an MR gradient waveform subject to a number of constraints.
% In particular it maximizes the b-value of a diffusion encoding
% pulse sequence subject to constraints on power, maximum gradient and maximum slew rate.
%
% If you use this in your research, please cite the following paper:
% Jens Sjölund, Filip Szczepankiewicz, Markus Nilsson, Daniel Topgaard, Carl-Fredrik Westin, Hans Knutsson,
% "Constrained optimization of gradient waveforms for generalized diffusion encoding",
% Journal of Magnetic Resonance, Volume 261, December 2015, Pages 157-168, ISSN 1090-7807,
% http://dx.doi.org/10.1016/j.jmr.2015.10.012.
% (http://www.sciencedirect.com/science/article/pii/S1090780715002451)
%
% If you use asymmetric waveforms with Maxwell compensation, please cite the following abstract (or later paper):
% Filip Szczepankiewicz and Markus Nilsson
% "Maxwell-compensated waveform design for asymmetric diffusion encoding"
% ISMRM 2018, Paris, France
% Download PDF at: https://goo.gl/vVGQq2
%
%
% Written by Jens Sjölund, jens.sjolund@elekta.com
% Maxwell compensation by Filip Szczepankiewicz, filip.szczepankiewicz@med.lu.se

%% Initialize parameters
if nargin == 0
    problem = optimizationProblem();
elseif nargin == 1
    problem = optimizationProblem(varargin{1});
else
    error('Invalid number of input arguments')
end


%% Set optimization parameters
options = optimoptions('fmincon','Algorithm','sqp',...
    'DerivativeCheck','off','Display','off',...
    'GradObj','on','GradConstr','on','MaxFunEval',1e5,'MaxIter',5e3);
warning('off', 'optimlib:fmincon:ConvertingToFull'); %Disables warning when SQP converts sparse matrices to full

%% Set up constraints
[A, b, firstDerivativeMatrix, secondDerivativeMatrix] = defineLinearInequalityConstraints(problem.N, problem.gMaxConstraint, problem.sMaxConstraint, problem.useMaxNorm);

[Aeq, beq] = defineLinearEqualityConstraints(problem.N, problem.zeroGradientAtIndex, problem.enforceSymmetry, firstDerivativeMatrix);

% Define nonlinear inequality constraints
nonlconFileName = getNonLinearConstraintsFileName(problem.N, problem.useMaxNorm);
if ~exist(nonlconFileName,'file')
    createConstraintGradientFunction(problem.N,problem.useMaxNorm); %Uses the symbolic toolbox to derive Jacobian ,SLOW!
end

%% Optimize
optimizationSuccess = false;
iter = 1;
while ~optimizationSuccess && iter <= 10
    
    x0 = getInitialGuess(problem, iter);
    dispInfo(problem, iter)
    
    tic
    
	[x,fval,exitflag,output,lambda,grad]  = fmincon(@(x) objFun(x), x0, A,b,Aeq,beq,[],[],@(x) feval(nonlconFileName,x,problem.tolIsotropy, ...
											problem.gMaxConstraint, problem.integralConstraint,problem.targetTensor, problem.tolMaxwell, ...
											problem.signs),options);
	
    optimizationTime = toc;
    
    disp(['Optimization took ' num2str(optimizationTime, 3) 's.']);
    
    optimizationSuccess = (exitflag > 0);
    
    iter = iter +1;
end

%% Evaluate and store results
gamma = 42.6e6*2*pi;
q = reshape(x(1:3*problem.N),[problem.N,3]);
g = [zeros(1,3);firstDerivativeMatrix*reshape(q,[problem.N 3]);zeros(1,3)]/problem.dt;
slew = [zeros(1,3);secondDerivativeMatrix*reshape(q,[problem.N 3]);zeros(1,3)]/(problem.dt)^2;
q = gamma*1e-6*q; %SI-units
q0 = reshape(x0(1:3*problem.N),[problem.N,3]);
q0 = gamma*1e-6*q0;
B = problem.dt*1e-3*(q'*q);
b = trace(B)*1e-6;%s/mm^2
C = b/(1e-6*gamma^2*(problem.gMaxConstraint*1e-3/problem.dt)^2*(problem.totalTimeActual*1e-3)^3);
kappa = 4*C;

etaOpt = problem.dt*max(diag(g'*g))/(problem.gMaxConstraint^2*problem.totalTimeActual);

result.q = q;
result.g = g;
result.slew = slew;
result.b = b;
result.B = B;
result.q0 = q0;
result.kappa = kappa;
result.etaOpt = etaOpt;
result.optimizerOutput = output;

if problem.doMaxwellComp
    rf = problem.signs;                
else
    rf = ones(size(problem.signs));
end

result.zind = problem.zeroGradientAtIndex;       % Keep this info for save function
result.rf   = [rf(1); rf; rf(end)];              % Spin direction
result.gwf  = bsxfun(@times, result.rf, g/1000); % T/m
result.dt   = problem.dt/1000;                   % s

end

function [A, b, firstDerivativeMatrix, secondDerivativeMatrix] = defineLinearInequalityConstraints(N, gMaxConstraint, sMaxConstraint, useMaxNorm)
firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1); % Center difference, shifted forward by half a step. Ghost points implemented as zero rows.
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
firstDerivativeMatrix = sparse(firstDerivativeMatrix); %SQP doesn't take advantage of this

if useMaxNorm == true %This is if we want to use max-norm on the gradients
    A1 = kron(eye(3),firstDerivativeMatrix);
    A1 = [A1 zeros(size(A1,1),1)]; %Add column of zeros for s
    b1 = gMaxConstraint*ones(size(A1,1),1);
else
    A1 = [];
    b1 = [];
end

%Constraint on change in gradients (slew rate)
secondDerivativeMatrix = diag(ones(N-1,1),-1)-2*diag(ones(N,1))+diag(ones(N-1,1),1);
secondDerivativeMatrix = sparse(secondDerivativeMatrix);%SQP doesnt take advantage of this
A2 = kron(eye(3),secondDerivativeMatrix);
A2 = [A2 zeros(size(A2,1),1)]; %Add column of zeros for s
b2 = sMaxConstraint*ones(size(A2,1),1);

A = [A1;-A1;A2;-A2]; %abs(Ax)<=b is equivalent to -Ax<=b && Ax<=b
b = [b1;b1;b2;b2];
end

function [Aeq, beq] = defineLinearEqualityConstraints(N, zeroGradientAtIndex, enforceSymmetry, firstDerivativeMatrix)
Aeq = zeros(2+length(zeroGradientAtIndex),N);
% Require start and end in q-space origin (echo condition)
Aeq(1,1)=1;
Aeq(2,N)=1;
% Require zero gradient at the specified indices
Aeq(2+(1:length(zeroGradientAtIndex)),:) = firstDerivativeMatrix(zeroGradientAtIndex,:);

% Enforce symmetry about zero gradient interval
if enforceSymmetry == true
    if isempty(zeroGradientAtIndex)
        indicesBefore = 1:floor(N/2);
        indicesAfter = floor(N/2) + (1:ceil(N/2));
    else
        indicesBefore = 1:zeroGradientAtIndex(1);
        indicesAfter = (zeroGradientAtIndex(end)+1):N;
    end
    
    assert(length(indicesBefore) == length(indicesAfter),'Cannot enforce symmetry since the number of time samples before and after zero gradient interval is not equal.')
    
    Nactivated = length(indicesBefore);
    Aeq = [Aeq;fliplr(eye(Nactivated)), zeros(Nactivated,N-2*Nactivated),-eye(Nactivated)]; %q(1) = q(end) and so on.
end

Aeq = kron(eye(3),Aeq);
Aeq = sparse([Aeq zeros(size(Aeq,1),1)]); %Add column of zeros for s
beq = zeros(size(Aeq,1),1);
end

function x0 = getInitialGuess(optimization, iter)
if strcmp(optimization.initialGuess,'user-provided')
    if iter > 1
        warning('Optimization with user-provided initial guess seems to have failed. Trying with random initialization instead.')
        x0 = randn(3*optimization.N+1,1);
    else
        x0 = optimization.x0; % Return user-specified initial guess.
    end
else
    if ~strcmp(optimization.initialGuess,'random') && iter == 1
        warning('Unrecognized initial guess string. Defaulting to random initialization.')
    end
    x0 = randn(3*optimization.N+1,1);
end
end

function dispInfo(problem, iter)
if iter == 1
    disp(['Optimizing ' problem.name])
else
    disp(['Optimizing ' problem.name ', attempt ' num2str(iter)])
end
end

