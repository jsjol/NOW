function [result, problem] = NOW_RUN(problem)
% function [result, problem] = NOW_RUN(problem)

% This function calls the private function "optimize". The repackaging is
% done so that "NOW_RUN" can be acessed by calls from anywhere in the
% matlab structure (unlike optimize, which is restricted).
[result, problem] = optimize(problem);