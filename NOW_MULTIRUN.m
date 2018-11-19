function [result, problem] = NOW_MULTIRUN(problem, iter)
% function [result, problem] = NOW_MULTIRUN(problem, iter)
% Perform multiple optimizations and keep the one that yields the highest
% b-value.

% This function calls the private function "optimize". The repackaging is
% done so that "NOW_RUN" can be acessed by calls from anywhere in the
% matlab structure (unlike optimize, which is restricted).


b_best = 0;

for i = 1:iter
    
    fprintf(['Running iteration ' num2str(i) '/' num2str(iter) '...'])
    
    [r, p] = NOW_RUN(problem);
    
    b_curr = r.b;
    
    if b_best < b_curr
        b_best = b_curr;
        result = r;
        problem = p;
    end
    
    fprintf(['Done! (b = ' num2str(b_curr, 2) ')\n'])
    
end