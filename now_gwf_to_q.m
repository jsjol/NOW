function q = now_gwf_to_q(gwf, dt)
% function q = now_gwf_to_q(gwf, dt)
% 
% Assume gwf is given as the effective waveform.

q = 2 * pi * now_gamma * cumsum(gwf, 1) * dt;

