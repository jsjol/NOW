function createConstraintGradientFunction(N,useMaxNorm)
q = sym('q',[3*N,1]);
targetTensor = sym('targetTensor',[3,3]);
signs = sym('signs', [N-1,1], 'real');

syms s B gMax tolIsotropy integralConstraint c1 c2 tolMaxwell real
integrationMatrix = eye(N);
integrationMatrix(1,1) = 0.5;
integrationMatrix(N,N) = 0.5;

B = (reshape(q,[N 3])'*integrationMatrix*reshape(q,[N 3]))*1/(N-1);
c1 = norm(B-s*targetTensor,'fro')^2-(s*tolIsotropy)^2;

firstDerivativeMatrix = -diag(ones(N,1))+diag(ones(N-1,1),1); % Center difference, shifted forward by half a step
firstDerivativeMatrix = firstDerivativeMatrix(1:end-1,:);
firstDerivativeMatrix = sparse(firstDerivativeMatrix);
g = firstDerivativeMatrix*reshape(q,[N 3]); %No need to include the zeros at start and end

%Power constraint: constrains the integral of g(t)^2. Since the gradients
%already represent the average of each time interval the trapezoid rule
%should not be used
c3 = g(:,1)'*g(:,1)-integralConstraint;
c4 = g(:,2)'*g(:,2)-integralConstraint;
c5 = g(:,3)'*g(:,3)-integralConstraint;

% Constraint for compensation of Maxwell terms (concomitant fields)
% For theory, please read/cite the following abstract (or later paper):
% Szczepankiewicz F, Westin, C?F, Nilsson M. Maxwell?compensated design 
% of asymmetric gradient waveforms for tensor?valued diffusion encoding. 
% Magn Reson Med. 2019;00:1–14. https://doi.org/10.1002/mrm.27828

M = g'*diag(signs)*g; % In MATLAB >= 2016b, broadcasting is preferable.
c6 = norm(M, 'fro') - tolMaxwell;


if useMaxNorm == false
    c2 = (sum(g.^2,2)-gMax^2)';%Nonlinear inequality constraint <= 0
    
    c = [c1 c2 c3 c4 c5 c6];
    gradc = jacobian(c,[q;s]).'; % transpose to put in correct form
else
    c = [c1 c3 c4 c5 c6];
    gradc = jacobian(c,[q;s]).'; % transpose to put in correct form
end

if useMaxNorm
    fileName = ['private/nonlcon' num2str(N) 'pointsMaxNorm'];
else
    fileName = ['private/nonlcon' num2str(N) 'points2Norm'];
end
matlabFunction(c,[],gradc,[],'file',fileName,'vars',{[q;s],tolIsotropy,gMax, integralConstraint,targetTensor, tolMaxwell, signs}); %The [] are for the inequality constraints


