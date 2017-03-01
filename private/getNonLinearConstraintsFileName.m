function nonlconFileName = getNonLinearConstraintsFileName(N, useMaxNorm)
if useMaxNorm == true
    nonlconFileName = ['nonlcon' num2str(N) 'pointsMaxNorm'];
else
    nonlconFileName = ['nonlcon' num2str(N) 'points2Norm'];
end
end