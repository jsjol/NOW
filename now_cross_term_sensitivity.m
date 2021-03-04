function [n, nt, H] = now_cross_term_sensitivity(result)
% function [n, nt, H] = now_cross_term_sensitivity(result)

gwf = result.g;
rf = result.rf;
dt = result.dt;

qt = now_gwf_to_q(gwf, dt);

H = now_gamma * cumsum(rf)*dt;
nt = cumsum(qt.*H)*dt;
n = nt(end,:);
