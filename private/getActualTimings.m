function [durationFirstPartActual, durationZeroGradientActual, durationSecondPartActual, totalTimeActual, zeroGradientAtIndex] = ...
    getActualTimings(durationFirstPartRequested, durationZeroGradientRequested, durationSecondPartRequested, discretizationSteps)
totalTime = durationFirstPartRequested + durationSecondPartRequested + durationZeroGradientRequested;
dt = totalTime/discretizationSteps;

startZeroGradientsIndex  = floor(durationFirstPartRequested/dt);
startSecondPartIndex   = floor((durationFirstPartRequested+durationZeroGradientRequested)/dt);

durationFirstPartActual = startZeroGradientsIndex * dt;
durationSecondPartActual = (discretizationSteps-startSecondPartIndex)* dt;
durationZeroGradientActual  = (startSecondPartIndex-startZeroGradientsIndex) * dt;
totalTimeActual = durationFirstPartActual + durationSecondPartActual + durationZeroGradientActual;

if durationZeroGradientActual > 0
    zeroGradientAtIndex = (startZeroGradientsIndex:startSecondPartIndex);
else
    zeroGradientAtIndex = [];
end

end

