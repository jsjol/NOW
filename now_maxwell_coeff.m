function [k_matrix, Maxwell_index, g2t] = now_maxwell_coeff(optimization)
% function [k_matrix, Maxwell_index, g2t] = now_maxwell_coeff(optimization)
%
% Calculation of Maxwell term and Maxwell index
%
% If you use asymmetric waveforms with Maxwell compensation, please
% cite the following abstract (or later paper):
% Filip Szczepankiewicz and Markus Nilsson
% "Maxwell-compensated waveform design for asymmetric diffusion encoding"
% ISMRM 2018, Paris, France
% Download PDF at: https://goo.gl/vVGQq2
%
% For theoretical background, see:
% Baron et al., The effect of concomitant gradient fields on diffusion 
% tensor imaging. Magn Reson Med, 2012. 68(4): p. 1190-201.

rf = optimization.result.rf;        % Spin direction
g  = optimization.result.g*1e-3;  	% T/m
dt = optimization.problem.dt*1e-3; 	% s
B0 = 1;                             % T. 
% Preferably, B0 would be the actual main magnetic field, but if 
% the actual value is higher than B0, we get a conservative estimate. 

M = (g'*diag(rf)*g)*dt; % In MATLAB >= 2016b, broadcasting is preferable.
Maxwell_index = norm(M, 'fro');

g2t = get_g2t(g', rf');

k_matrix = get_k_matrix(g, rf, B0, dt);

end

function g2t = get_g2t(g, rf)
    gxx = cumsum(g(1,:).*g(1,:).*rf);
    gyy = cumsum(g(2,:).*g(2,:).*rf);
    gzz = cumsum(g(3,:).*g(3,:).*rf);
    gxy = cumsum(g(1,:).*g(2,:).*rf);
    gxz = cumsum(g(1,:).*g(3,:).*rf);
    gyz = cumsum(g(2,:).*g(3,:).*rf);

    g2t = [gxx; gyy; gzz; gxy; gxz; gyz];
end

function k = get_k_matrix(g, rf, B0, dt)
gammaOver2Pi = now_gamma; % 1 / (T s)

% Waveform moment vector scale matrix
% k-vector along direction [x y z] for one part of the waveform (pre or
% post) is
%
%              |  GzGz    0      -2GxGz        |
% k = integral(|  0       GzGz   -2GyGz        |) * [x y z]' / (4B0) * gamma * dt
%              | -2GxGz  -2GyGz   4(GxGx+GyGy) |
%
% Final moment is k = k_pre - k_post

gx = g(:, 1);
gy = g(:, 2);
gz = g(:, 3);

v = [gz'*(rf.*gz), 0, -2*gx'*(rf.*gz);
     0, gz'*(rf.*gz), -2*gy'*(rf.*gz);
     -2*gx'*(rf.*gz), -2*gy'*(rf.*gz), 4*(gx'*(rf.*gx) + gy'*(rf.*gy))];
 
k = gammaOver2Pi * v / (4 * B0) * dt;

end
