function [c,ceq,gradc,gradceq] = nonlconAnalytic(x, problem)
%NONLCONANALYTIC evaluates the nonlinear constraints and their gradients.

q = x(1:end-1);
s = x(end);
N = problem.N;
Q = reshape(q,[N 3]);

% Upsample
% upsampled_N = 500;
% interpolation_matrix = interp1(1:N, eye(N), linspace(1, N, upsampled_N), 'cubic');
% Q = interpolation_matrix * Q;%interp1(1:N, Q, 1:0.5:N);
% N = size(Q, 1);
% problem.signs = interp1(1:length(problem.signs), problem.signs, linspace(1, length(problem.signs), upsampled_N - 1)', 'nearest');

% Remove action on axes that are turned off - disabled awaiting further
% tests or improved handling of disabled axes
Q = Q .* ((diag(problem.targetTensor))>0)';

integrationWeights = ones(N,1);
integrationWeights(1) = 0.5;
integrationWeights(N) = 0.5;
integrationWeights = integrationWeights / (N-1);
weightedQ = bsxfun(@times, Q, integrationWeights);
B = Q'*weightedQ;

% First constraint: tensor encoding,
% the first term is the squared Frobenius norm of the deviation
c1 = trace((B-s*problem.targetTensor)'*(B-s*problem.targetTensor))-(s*problem.tolIsotropy)^2;

% Gradient of tensor encoding constraint
dc1_dB = 2*B - 2*s*problem.targetTensor;
% B is a "matrix quadratic form", with a derivative that has two
% identical but transposed terms
firstTerm = kron(eye(3), weightedQ');
secondTerm = reshape(firstTerm, [3,3, 3*N]);
secondTerm = permute(secondTerm, [2, 1, 3]);
secondTerm = reshape(secondTerm, [9, 3*N]);
dB_dq = firstTerm + secondTerm;

dc1_dq = reshape(dc1_dB, [1, 9]) * dB_dq;
dc1_ds = 2*s*trace(problem.targetTensor'*problem.targetTensor) - 2*trace(B'*problem.targetTensor) - 2*s*problem.tolIsotropy^2;
dc1_dx = [dc1_dq, dc1_ds];

firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1); % Center difference, shifted forward by half a step
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
firstDerivativeMatrix = sparse(firstDerivativeMatrix);
firstDerivativeMatrix = firstDerivativeMatrix / problem.dt;
g = firstDerivativeMatrix*Q; %No need to include the zeros at start and end

% Constrain optimization norm
if problem.useMaxNorm
    c2 = [];
    dc2_dx = [];
    
else % L2-norm
    c2 = (sum(g.^2,2)-problem.gMax^2)'; %Nonlinear inequality constraint <= 0
    
%     dc2_dq1 = [bsxfun(@times, 2*firstDerivativeMatrix, g(:, 1))...
%         bsxfun(@times, 2*firstDerivativeMatrix, g(:, 2))...
%         bsxfun(@times, 2*firstDerivativeMatrix, g(:, 3))];
    
    % Steven Baete, below is faster (~x4)
    dc2_dq = 2*[firstDerivativeMatrix.*repmat(g(:,1),[1,N]),...
            firstDerivativeMatrix.*repmat(g(:,2),[1,N]),...
            firstDerivativeMatrix.*repmat(g(:,3),[1,N])];
      
    dc2_dx = [dc2_dq, zeros(N-1, 1)];
    
end

%Power constraint: constrains the integral of g(t)^2. Since the gradients
%already represent the average of each time interval the trapezoid rule
%should not be used
c3 = g(:,1)'*g(:,1)-problem.integralConstraint;
c4 = g(:,2)'*g(:,2)-problem.integralConstraint;
c5 = g(:,3)'*g(:,3)-problem.integralConstraint;

dc3_dx = zeros(1, 3*N+1);
dc3_dx(1:N) = 2 * g(:,1)' * firstDerivativeMatrix;

dc4_dx = zeros(1, 3*N+1);
dc4_dx((N+1):(2*N)) = 2 * g(:,2)' * firstDerivativeMatrix;

dc5_dx = zeros(1, 3*N+1);
dc5_dx((2*N+1):(3*N)) = 2 * g(:,3)' * firstDerivativeMatrix;


% Constraint for compensation of Maxwell terms (concomitant fields)
% For theory, please read/cite the following paper:
% Szczepankiewicz F, Westin, C-F, Nilsson M. Maxwell-compensated design
% of asymmetric gradient waveforms for tensor-valued diffusion encoding.
% Magn Reson Med. 2019, Vol. 82, Issue 4, p.1424-1437.
% https://doi.org/10.1002/mrm.27828

if isinf(problem.MaxwellIndex)
    c6 = [];
    dc6_dx = [];
else
    signedg = bsxfun(@times, g, problem.signs);
    M = g'*signedg;
    m = sqrt(trace(M'*M)); % 'Maxwell index'
    c6 = m - problem.MaxwellIndex; % Check whether to square or not (as in ref)
    
    dc6_dM = 1/m * M;
    % M is a "matrix quadratic form", so we use the same procedure as above,
    % but for performance reasons we keep it inline instead of
    % defining a function
    weightedQ = firstDerivativeMatrix'*signedg;
    firstTerm = kron(eye(3), weightedQ');
    secondTerm = reshape(firstTerm, [3,3, 3*N]);
    secondTerm = permute(secondTerm, [2, 1, 3]);
    secondTerm = reshape(secondTerm, [9, 3*N]);
    dM_dq = firstTerm + secondTerm;
    dc6_dq = reshape(dc6_dM, [1, 9]) * dM_dq;
    dc6_dx = [dc6_dq, 0];
end


% Constraint for ballistic motion encoding
% For theory, please read/cite the following paper:
% Szczepankiewicz F, Sjolund J, ... Westin C-F.
% Motion-compensated gradient waveforms for tensor-valued
% diffusion encoding by constrained numerical optimization
% Magn Reson Med. 2020
% https://doi.org/10.1002/mrm.28551

% Non-linear constraints for motion compensation
nonlinear_ind = find(~problem.motionCompensation.linear);
if isempty(nonlinear_ind)
    c7 = [];
    dc7_dx = [];
else
    for i = 1:length(nonlinear_ind)
        c7 = zeros(1, length(nonlinear_ind));
        dc7_dx = zeros(length(nonlinear_ind), 3*N + 1);
        
        % dt is in ms; % Behaves better if calculation uses ms but
        % requested tolerance has some strange units. Fixed this by
        % rescaling the maxMagnitude by 1000^order.
        
        t = ((1:N)-1/2) * problem.dt;
        
        gamma = 2.6751e+08; % radians / T / s for hydrogen.
        
        order = problem.motionCompensation.order(nonlinear_ind(i));
        moment_weighting = - order * problem.dt * t.^(order-1);
        moment_vector = moment_weighting * Q;
        c7(i) = sum(moment_vector.^2) - (problem.motionCompensation.maxMagnitude(nonlinear_ind(i)) * 1000^order / (gamma * 1e-6))^2;
        dc7_dx(i, 1:(3*N)) = 2 * kron(moment_vector, moment_weighting);
    end
end


ceq     = [];
gradceq = [];
c       = [c1 c2 c3 c4 c5 c6 c7];
gradc   = [
    dc1_dx;
    dc2_dx;
    dc3_dx;
    dc4_dx;
    dc5_dx;
    dc6_dx;
    dc7_dx]'; % transpose the Jacobian to put in fmincon form

% Downsample to original time resolution
% gradc_s = gradc(end, :);
% gradc = [kron(eye(3), interpolation_matrix') * gradc(1:end-1, :); 
%          gradc_s];
