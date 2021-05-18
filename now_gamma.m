function gamma = now_gamma(nuc)
% function gamma = now_gamma(nuc)
%
% Returns the gyromagnetic constant for given nucleus in units of 1/s/T.

if nargin < 1
    nuc = 'H';
end


switch nuc
    case 'H'
        gamma = 42.6e6;
        
    otherwise
        error('Nucleus not recognized')
end

