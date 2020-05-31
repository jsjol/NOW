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
% If you use asymmetric waveforms with Maxwell compensation, 
% please cite the following abstract (or later paper):
% Szczepankiewicz F, Westin, C-F, Nilsson M. Maxwell-compensated design 
% of asymmetric gradient waveforms for tensor-valued diffusion encoding. 
% Magn Reson Med. 2019;00:1-14. https://doi.org/10.1002/mrm.27828
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
    'DerivativeCheck','off','FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-4,...% Gradient check (barely) fails, likely for numerical reasons because errors increase with smaller step sizes than 1e-5.
    'Display','off', 'GradObj','on','GradConstr','on','MaxFunEval', ...
    problem.MaxFunEval, 'MaxIter', problem.MaxIter);
warning('off', 'optimlib:fmincon:ConvertingToFull'); %Disables warning when SQP converts sparse matrices to full

%% Set up constraints
[A, b] = defineLinearInequalityConstraints(problem);

[Aeq, beq] = defineLinearEqualityConstraints(problem);

% Define nonlinear inequality constraints

% Uncomment when using symbolic graidents
% nonlconFileName = getNonLinearConstraintsFileName(problem.N, problem.useMaxNorm);
% if ~exist(nonlconFileName,'file')
%     createConstraintGradientFunction(problem.N,problem.useMaxNorm); %Uses the symbolic toolbox to derive Jacobian ,SLOW!
% end

%% Optimize
optimizationSuccess = false;
iter = 1;
while ~optimizationSuccess && iter <= 10
    
    x0 = getInitialGuess(problem, iter);
    dispInfo(problem, iter)
    
    tic
    % Symbolic gradients
% 	[x,fval,exitflag,output,lambda,grad]  = fmincon(@(x) objFun(x), x0, A,b,Aeq,beq,[],[],@(x) feval(nonlconFileName,x,problem.tolIsotropy, ...
% 											problem.gMaxConstraint, problem.integralConstraint,problem.targetTensor, problem.tolMaxwell*problem.dt^2, ...
% 											problem.signs),options);
    % Analytic gradients                                  
    [x,fval,exitflag,output,lambda,grad]  = fmincon(@(x) objFun(x), x0, A,b,Aeq,beq,[],[], @(x) nonlconAnalytic(x,problem.tolIsotropy, ...
											problem.gMaxConstraint, problem.integralConstraint,problem.targetTensor, problem.tolMaxwell*problem.dt^2, ...
											problem.signs, problem.useMaxNorm, problem.motionCompensation, problem.dt),options);
	
    optimizationTime = toc;
    
    disp(['Optimization took ' num2str(optimizationTime, 3) 's.']);
    
    optimizationSuccess = (exitflag > 0);
    
    iter = iter +1;
    
    if ~problem.redoIfFailed
        if ~optimizationSuccess
            disp('Optimization failed but will not be repeated!')
        end
        break
    end
end

%% Evaluate and store results
[firstDerivativeMatrix, secondDerivativeMatrix] = getDerivativeMatrices(problem);

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
result.optimizationTime = optimizationTime;
result.iter = iter-1;
result.rawq = x(1:(end-1));

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

function [firstDerivativeMatrix, secondDerivativeMatrix] = getDerivativeMatrices(problem)
    firstDerivativeMatrix = -diag(ones(problem.N,1))+diag(ones(problem.N-1,1),1); % Center difference, shifted forward by half a step. Ghost points implemented as zero rows.
    firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
    firstDerivativeMatrix = sparse(firstDerivativeMatrix); %SQP doesn't take advantage of this

    secondDerivativeMatrix = diag(ones(problem.N-1,1),-1)-2*diag(ones(problem.N,1))+diag(ones(problem.N-1,1),1);
    secondDerivativeMatrix = sparse(secondDerivativeMatrix);%SQP doesnt take advantage of this
end

function [A, b] = defineLinearInequalityConstraints(problem)

    [firstDerivativeMatrix, secondDerivativeMatrix] = getDerivativeMatrices(problem);

    if problem.useMaxNorm == true %This is if we want to use max-norm on the gradients
        A1 = kron(eye(3),firstDerivativeMatrix);
        A1 = [A1 zeros(size(A1,1),1)]; %Add column of zeros for s
        b1 = problem.gMaxConstraint*ones(size(A1,1),1);
    else
        A1 = [];
        b1 = [];
    end

    %Constraint on change in gradients (slew rate)

    A2 = kron(eye(3),secondDerivativeMatrix);
    A2 = [A2 zeros(size(A2,1),1)]; %Add column of zeros for s
    b2 = problem.sMaxConstraint*ones(size(A2,1),1);

    % Motion compensation - linear formulation
%     motion_moment = zeros(length(problem.motionCompensation.order), 3 * problem.N);
%     for i = 1:length(problem.motionCompensation.order)
%        n = problem.motionCompensation.order(i);
%        t = ((1:problem.N)-1/2) * problem.dt;
%        n_moment = - n * t.^(n-1) * problem.dt; % Single channel
%        %motion_moment(i, :) = % max norm?
%     end
%     motion_moment = [motion_moment, zeros(size(motion_moment,1), 1)]; %Add column of zeros for s

    A = [A1;-A1;A2;-A2]; %abs(Ax)<=b is equivalent to -Ax<=b && Ax<=b
    b = [b1;b1;b2;b2];
end

function [Aeq, beq] = defineLinearEqualityConstraints(problem)

[firstDerivativeMatrix, ~] = getDerivativeMatrices(problem);

Aeq = zeros(2+length(problem.zeroGradientAtIndex),problem.N);
% Require start and end in q-space origin (echo condition)
Aeq(1,1)=1;
Aeq(2,problem.N)=1;
% Require zero gradient at the specified indices
Aeq(2+(1:length(problem.zeroGradientAtIndex)),:) = firstDerivativeMatrix(problem.zeroGradientAtIndex,:);

% Enforce symmetry about zero gradient interval
if problem.enforceSymmetry == true
    if isempty(problem.zeroGradientAtIndex)
        indicesBefore = 1:floor(problem.N/2);
        indicesAfter = floor(problem.N/2) + (1:ceil(problem.N/2));
    else
        indicesBefore = 1:problem.zeroGradientAtIndex(1);
        indicesAfter = (problem.zeroGradientAtIndex(end)+1):problem.N;
    end
    
    assert(length(indicesBefore) == length(indicesAfter),'Cannot enforce symmetry since the number of time samples before and after zero gradient interval is not equal.')
    
    Nactivated = length(indicesBefore);
    Aeq = [Aeq;fliplr(eye(Nactivated)), zeros(Nactivated,problem.N-2*Nactivated),-eye(Nactivated)]; %q(1) = q(end) and so on.
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

