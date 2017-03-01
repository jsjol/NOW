function result = scripted_NOW_Example()
%Example of scripted call to Numerical Optimization of gradient Waveform
%(NOW). 
%
%MATLAB requires it to be a function in order to access the
%functions in the private folder.

%Change the parameters below to your liking. Those not
%specified are default-initialized as follows:
% Max gradient = 80 milliTesla/m
% Max slew rate = 100 milliTesla/m/milliSecond = 100 T/m/s
% Pulse-time = 50 milliSecond
% Eta (heat dissipation parameter) = 1
% Discretization points = 50
% Target tensor = eye(3)
% Initialguess = 'random'
% enforceSymmetry = false;
% redoIfFailed = true;
% name = 'NOW'

% Written by Jens Sjölund and Filip Szczepankiewicz

%%  PREP

requestedGradientAmplitude = 80; % In mT/m
requestedSlewRate = 60; % In T/(sm)

% Request encoding and pause times based on sequence [ms]
settings.durationFirstPartRequested = 51;
settings.durationSecondPartRequested = 40;
settings.durationZeroGradientRequested = 8;

settings.targetTensor = [1 0 0;0 1 0;0 0 1];
settings.N = 100;
settings.gMax = requestedGradientAmplitude;
settings.sMax = requestedSlewRate;
settings.eta = 0.5; %In range (0,1]
settings.initialGuess = 'random';
settings.redoIfFailed = true;
settings.enforceSymmetry = false;

problem = optimizationProblem(settings); 

%% PRINT REQUESTED AND TRUE TIMES
clc
fprintf(1, '--------- Requested timing parameters: --------- \n');
fprintf(1, 'Dur1 = %7.3f  Dur2 = %7.3f  DurPi = %6.3f  [ms]\n\n', problem.durationFirstPartRequested, problem.durationSecondPartRequested, problem.durationZeroGradientRequested);
fprintf(1, '---------   Actual timing parameters:  --------- \n');
fprintf(1, 'Dur1 = %7.3f  Dur2 = %7.3f  DurPi = %6.3f  [ms]\n', problem.durationFirstPartActual, problem.durationSecondPartActual, problem.durationZeroGradientActual);


%% RUN OPTIMIZATION

[result, problem] = optimize(problem); % You could also pass settings directly


%% PLOT RESULT
plot(0:problem.dt:problem.totalTimeActual,result.g)
xlabel('Time [ms]')
ylabel('Gradient amplitude [mT/m]')
measurementTensor = result.B


