function createConstraintGradientFunction(N,useMaxNorm)
q = sym('q',[3*N,1]);
targetTensor = sym('targetTensor',[3,3]);
s_vec = sym('s_vec', [N-1,1], 'real');

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

x = 1;
y = 2;
z = 3;

% Constraint for compensation of Maxwell terms (concomitant fields)
% For theory, please read/cite the following abstract (or later paper):
% Filip Szczepankiewicz and Markus Nilsson
% "Maxwell-compensated waveform design for asymmetric diffusion encoding"
% ISMRM 2018, Paris, France
% Download PDF at: https://goo.gl/vVGQq2

gxx = sum( g(:,x).*g(:,x).*s_vec );
gyy = sum( g(:,y).*g(:,y).*s_vec );
gzz = sum( g(:,z).*g(:,z).*s_vec );
gxy = sum( g(:,x).*g(:,y).*s_vec );
gxz = sum( g(:,x).*g(:,z).*s_vec );
gyz = sum( g(:,y).*g(:,z).*s_vec );

M = [gxx gxy gxz; gxy gyy gyz; gxz gyz gzz];

c6 = sqrt(trace(M*M)) - tolMaxwell;



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
matlabFunction(c,[],gradc,[],'file',fileName,'vars',{[q;s],tolIsotropy,gMax, integralConstraint,targetTensor, tolMaxwell, s_vec}); %The [] are for the inequality constraints


