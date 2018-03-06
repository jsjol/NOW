function [k, m, g2t] = now_maxwell_coeff(g, r)
% function [k, m, g2t] = now_maxwell_coeff(g, r)
%
% Quick and dirty implementation of Maxwell term calculation, and Maxwell index
%
% If you use asymmetric waveforms with Maxwell compensation, please cite the following abstract (or later paper):
% Filip Szczepankiewicz and Markus Nilsson
% "Maxwell-compensated waveform design for asymmetric diffusion encoding"
% ISMRM 2018, Paris, France
% Download PDF at: https://goo.gl/vVGQq2
%
% For theoretical background, see:
% Baron et al., The effect of concomitant gradient fields on diffusion 
% tensor imaging. Magn Reson Med, 2012. 68(4): p. 1190-201.

g  = g'*1e-3;               % T/m
dt = r.problem.dt*1e-3;     % s
B0 = 1;                     % T, worst case scenario (?)
zi = r.result.optimizerProblem.zeroGradientAtIndex;
N  = r.problem.N;

if isempty(zi)
    ga = g;
    gb = zeros(3,1);
    rf = ones(1, size([ga gb], 2));
    
else
    coi = zi(ceil(numel(zi)/2));
    
    pre  = 1:coi;
    post = (coi+1):(N+1);
    
    ga = g(:, pre);
    gb = g(:, post);
    
    rf = ones(1, size(g, 2));
    rf(post:end) = -1;
    
end

gxx = cumsum(g(1,:).*g(1,:).*rf);
gyy = cumsum(g(2,:).*g(2,:).*rf);
gzz = cumsum(g(3,:).*g(3,:).*rf);
gxy = cumsum(g(1,:).*g(2,:).*rf);
gxz = cumsum(g(1,:).*g(3,:).*rf);
gyz = cumsum(g(2,:).*g(3,:).*rf);

g2t = [gxx; gyy; gzz; gxy; gxz; gyz];

k = k_matrix(ga, gb, B0, dt);

m = maxwell_index(ga, gb, dt);

end


function m = maxwell_index(ga, gb, dt)

M = (ga*ga'-gb*gb')*dt;

m = sqrt(trace(M*M)); % (T/m)^2 s

end


function k = k_matrix(ga, gb, B0, dt)
gammaOver2PI = 42.6e6; % 1 / (T s)

% Waveform moment vector scale matrix
% k-vector along direction [x y z] for one part of the waveform (pre or
% post) is
%
%              |  GzGz    0      -2GxGz        |
% k = integral(|  0       GzGz   -2GyGz        |) * [x y z]' / (4B0) * gamma * dt
%              | -2GxGz  -2GyGz   4(GxGx+GyGy) |
%
% Final moment is k = k_pre - k_post

v = @(G)[...
    sum(G(3,:).^2,2),                                  0,                   -2*sum(G(1,:).*G(3,:),2)  ;
    0,                   sum(G(3,:).^2,2),                   -2*sum(G(2,:).*G(3,:),2)  ;
    -2*sum(G(1,:).*G(3,:),2),   -2*sum(G(2,:).*G(3,:),2),    4*(sum(G(1,:).^2,2) + sum(G(2,:).^2,2)) ];


k = gammaOver2PI * (v(ga) - v(gb)) / (4 * B0) * dt;
end
