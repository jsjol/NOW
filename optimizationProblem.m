classdef optimizationProblem
    %OPTIMIZATIONPROBLEM Store all properties needed to run NOW.
    %   Parameters not specified by the user are default-initialized as follows:
    %
    %   Max gradient = 80 milliTesla/m
    %   Max slew rate = 100 milliTesla/m/milliSecond = 100 T/m/s
    %   Pulse-time = 50 milliseconds
    %   Eta (heat dissipation parameter) = 1
    %   Discretization points = 50
    %   Target tensor = eye(3)
    %   Initialguess = 'random'
    %   zeroGradientAtIndex = [], i.e. only at start and end
    %   enforceSymmetry = false;
    %   redoIfFailed = true;
    
    properties (Access = public)
        targetTensor = 1/3*eye(3); %Isotropic encoding tensor
        N = 50;
        initialGuess = 'random';
        useMaxNorm = true;
        gMax = 80;
        sMax = 100;
        durationFirstPartRequested = 25;
        durationSecondPartRequested = 25;
        durationZeroGradientRequested = 0;
        eta = 1;
        enforceSymmetry = false;
        redoIfFailed = true;
        name = 'NOW';
    end
    
    properties (SetAccess = private)
        zeroGradientAtIndex = [];
        tolIsotropy = 1e-4;
        durationFirstPartActual
        durationZeroGradientActual
        durationSecondPartActual
        totalTimeActual
        dt
        gMaxConstraint
        sMaxConstraint
        integralConstraint
    end
    
    methods (Access = public)
        function obj = optimizationProblem(varargin)
            if nargin > 0
                settings = varargin{1};
                
%                 if isa(settings, 'optimizationProblem')
%                     obj = settings;
%                     return
%                 end
                
                % Overwrite defaults with user-specified settings
                fieldNames = fieldnames(settings);
                for i = 1:length(fieldNames)
                    eval(['obj.' fieldNames{i} ' = getfield(settings, fieldNames{i});'])
                end
            end
            
            % Get actual times after discretization
            [obj.durationFirstPartActual, obj.durationZeroGradientActual, obj.durationSecondPartActual, obj.totalTimeActual, obj.zeroGradientAtIndex] = ...
                    getActualTimings(obj.durationFirstPartRequested, obj.durationZeroGradientRequested, obj.durationSecondPartRequested, obj.N);
            
            
            % Compute private variables
            obj.dt = obj.totalTimeActual/obj.N; %Time step in milliseconds. Division by N instead of N-1 due to half step shift in gradients.
            obj.gMaxConstraint = obj.gMax*obj.dt;
            obj.sMaxConstraint = obj.sMax*obj.dt^2; 
            obj.integralConstraint = obj.eta*obj.gMaxConstraint^2*obj.totalTimeActual/obj.dt;

        end
    end
    
end

