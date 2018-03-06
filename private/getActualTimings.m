function [durationFirstPartActual, durationZeroGradientActual, durationSecondPartActual, totalTimeActual, zeroGradientAtIndex] = ...
    getActualTimings(durationFirstPartRequested, durationZeroGradientRequested, durationSecondPartRequested, discretizationSteps, forceSymmetry)

if ~forceSymmetry
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
    
else
    
    durBothRequested = min([durationFirstPartRequested durationSecondPartRequested]);
    
    totalTime = 2*durBothRequested + durationZeroGradientRequested;
    dt = totalTime/discretizationSteps;
    
    num_zero =  ceil(durationZeroGradientRequested / dt);
    
    if isEven(discretizationSteps - num_zero)
        % OK, do nothing
    else
        num_zero = num_zero-1;
    end
    
    startZeroGradientsIndex  = (discretizationSteps - num_zero)/2;
    startSecondPartIndex     = (discretizationSteps - num_zero)/2 + num_zero;
    
    durationFirstPartActual = startZeroGradientsIndex * dt;
    durationSecondPartActual = (discretizationSteps-startSecondPartIndex)* dt;
    durationZeroGradientActual  = (startSecondPartIndex-startZeroGradientsIndex) * dt;
    totalTimeActual = durationFirstPartActual + durationSecondPartActual + durationZeroGradientActual;
    
    if durationFirstPartActual ~= durationSecondPartActual
        error('asfafasfa')
    end
    
    if durationZeroGradientActual > 0
        zeroGradientAtIndex = (startZeroGradientsIndex:startSecondPartIndex);
    else
        zeroGradientAtIndex = [];
    end
    
    
    
end


end

