classdef optimizationProblem
    %OPTIMIZATIONPROBLEM Store all properties needed to run NOW.
    %   Parameters not specified by the user are default-initialized as follows:
    %
    %   Max gradient = 80 milliTesla/m
    %   Max slew rate = 100 milliTesla/m/milliSecond = 100 T/m/s
    %   Eta (heat dissipation parameter) = 1
    %   Discretization points = 77
    %   Target tensor = eye(3)
    %   Initialguess = 'random'
    %   zeroGradientAtIndex = [], i.e. only at start and end
    %   enforceSymmetry = false;
    %   redoIfFailed = true;
    %   useMaxNorm = false;
    %   doMaxwellComp = true;
    %   MaxwellIndex = 100;
    %   Motion compensation: disabled (when enabled, magnitude unit is s^order / m)
    
    
    properties (Access = public)
        targetTensor = eye(3); % Isotropic encoding tensor
        N = 77;
        initialGuess = 'random';
        useMaxNorm = false;
        gMax = 80;
        sMax = 100;
        durationFirstPartRequested = 28;
        durationSecondPartRequested = 22;
        durationZeroGradientRequested = 8;
        eta = 1;
        enforceSymmetry = false;
        redoIfFailed = true;
        name = 'NOW';
        x0 = [];
        doMaxwellComp = true;
        MaxwellIndex = 100;
        MaxFunEval = 1e5;
        MaxIter    = 5e3;
        motionCompensation = struct('order', [], 'maxMagnitude', [], 'linear', [])
        doBackgroundCompensation = 0; % 0 = off; 1 = general timing cond.; 2 = specific timing cond.
        startTime = 0; % Time from the excitataion (t=0) to the first gradietn waveform sample in ms.
    end
    
    properties (SetAccess = private)
        zeroGradientAtIndex = [];
        tolIsotropy = .5e-2;
        tolMaxwell
        signs
        tolSlew
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
                
                % Overwrite defaults with user-specified settings
                fieldNames = fieldnames(settings);
                for i = 1:length(fieldNames)
                    eval(['obj.' fieldNames{i} ' = getfield(settings, fieldNames{i});'])
                end
            end
            
            % Get actual times after discretization
            [obj.durationFirstPartActual, obj.durationZeroGradientActual, obj.durationSecondPartActual, obj.totalTimeActual, obj.zeroGradientAtIndex] = ...
                getActualTimings(obj.durationFirstPartRequested, obj.durationZeroGradientRequested, obj.durationSecondPartRequested, obj.N, obj.enforceSymmetry);
            
            
            % Compute private variables
            obj.dt = obj.totalTimeActual/obj.N; %Time step in milliseconds. Division by N instead of N-1 due to half step shift in gradients.
            obj.gMaxConstraint = obj.gMax*obj.dt;
            obj.sMaxConstraint = obj.sMax*obj.dt^2;
            obj.integralConstraint = obj.eta*obj.gMaxConstraint^2*obj.totalTimeActual/obj.dt;
            obj.tolMaxwell = obj.MaxwellIndex/obj.dt;
            
            % Turn Maxwell compensation off if requested
            % Can also be written as a more confusing and compact form:
            % obj.tolMaxwell = obj.MaxwellIndex/obj.dt / (obj.doMaxwell>0);
            % But doMaxwell can also be entirely replaced by
            % obj.MaxwellIndex where < inf indicates that its ON but this will
            % complicate the GUI.
            if ~obj.doMaxwellComp
                obj.tolMaxwell = inf;
            end
            
            
            %% Create spin dephasing direction vector
            if ~isempty(obj.zeroGradientAtIndex)
                signs = ones(obj.N - 1,1); % Ghost points excluded during opt
                
                % Assume that sign change happens in the middle of the pause
                mi = median(obj.zeroGradientAtIndex);
                
                signs(round(mi):end) = -1;
                
                if (mi==round(mi))
                    signs(round(mi)) = 0;
                end
                
                obj.signs = signs;
            end
            
            
            %% Motion compensation
            if length(obj.motionCompensation.maxMagnitude) ~= length(obj.motionCompensation.order)
                error('motionCompensation.maxMagnitude must have the same size as motionCompensation.order.')
            end
            
            if isempty(obj.motionCompensation.maxMagnitude)
                obj.motionCompensation.linear = [];
            else
                % Infer empty motionCompensation.linear from values of
                % motionCompensation.maxMagnitude
                obj.motionCompensation.linear = (obj.motionCompensation.maxMagnitude <= 0);
            end
            
            %% Cross-term-compensation
            switch obj.doBackgroundCompensation
                case 0 % No compensation
                    % Do nothing
                    
                case 1 % General timing condition
                    % This case requires velocity compensation.
                    
                    ind_velo = find(obj.motionCompensation.order==1, 1);
                    
                    % Check if velocity compensation is requested at all.
                    % If not, force linear constraint.
                    % Else, check that requested velo compensation is
                    % linear. If not, throw error since instructions cannot
                    % be fulfilled.
                    
                    if isempty(ind_velo)
                        ind = numel(obj.motionCompensation.order)+1;
                        obj.motionCompensation.order(ind) = 1;
                        obj.motionCompensation.linear(ind) = 1;
                        obj.motionCompensation.maxMagnitude(ind) = 0;
                    else
                        if ~obj.motionCompensation.linear(ind_velo)
                            error('Cross-term-compensation for a general timing requires velocity compensation!')
                        end
                    end
                    
                case 2 % Specific timing condition
                    % Specific timing condition requires that the start
                    % time is set. The start time can be zero, but it is
                    % unlikely that this is an accurate setup, so we throw
                    % a warning.
                    
                    if obj.startTime == 0
                        warning('Start time for waveform is t = 0, which is an unlikely setting. Please check!')
                    end
                    
                    if obj.startTime < 0
                        error('Start time cannot be smaller than zero!')
                    end
                    
                otherwise
                    error('Selection for Cross-term-compensation not recognized! Use value 0, 1 or 2.')
            end
            
        end
    end
end

