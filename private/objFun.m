function [f ,g] = objFun(x)
f = -x(end); %fmincon minimizes the objective, so a minus sign is used to maximize instead
g = zeros(size(x));
g(end) = -1;

OptTime = toc;

if OptTime > (length(g)^2/.5e3)
    error(['Timed out after ' num2str(OptTime) ' seconds!']);
end