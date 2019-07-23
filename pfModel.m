% Place Field Model object
%     contains information needed for general place field model
%
% Properties:
%     rate - function that predicts firing rate from input variables
%            has the form rate(b, inputs) where B is an nParamsx1 array of
%            parameter values and INPUTS is a matrix where each column is
%            an input variable and each row is a time point (eg. [position,
%            phase, speed])
%     dt - length of each time bin, typically 1/sampling rate
%     nParams - number of parameters in model
%     fitBounds - upper and lower bounds on fitting in fmincon function
%     initBounds - range from which each parameter is randomly chosen with 
%                  uniform probability in multistart fitting
%
% Methods: Described in function defintions
%
% Author: Kathryn McClain
% Date: 7/23/2019
%

classdef pfModel
    
    properties
        nParams
        initBounds
        fitBounds
        rate
        dt
    end
    
    methods
        

        function obj = pfModel(rate, nParams, varargin)
            %%%% Instantiate model object %%%
            % see header
            
            %parseing inputs
            p = inputParser;
            addParameter(p,'initBounds',@isnumeric)
            addParameter(p,'fitBounds',@isnumeric)
            
            parse(p,varargin{:});
            
            obj.rate = rate;
            obj.nParams = nParams;
            obj.initBounds = p.Results.initBounds;
            obj.fitBounds = p.Results.fitBounds;
            
            if isempty(obj.initBounds)
                
                obj.initBounds = zeros(nParams,2);
                obj.initBounds(:,1) = -Inf;
                obj.initBounds(:,2) = Inf;
                
            end
            if isempty(obj.fitBounds)
                
                obj.fitBounds = zeros(nParams,2);
                obj.fitBounds(:,1) = -Inf;
                obj.fitBounds(:,2) = Inf;
            end
            
            obj.dt = 1/1250; %timescale, 1/sampling rate
        end
        

        function [param, LL] = fitModel(obj,data,initPoint)
            %%% Fit model to data %%%
            % estimates model parameters from place field data
            %
            % inputs:
            %   data - matrix of the form [timestamp, spike/no spike, input variables(1 or more)]
            %   initPoint - initial point in fitting, size nParamsx2
            %
            % outputs:
            %   param - array of parameters estimated, size nParamsx1
            %   LL - log-likelihood of estimated parameters producing spiking
            %       in data
            
            spike = data(:,2);
            inputVars = data(:,3:end);
            
            lowerBound = obj.fitBounds(:,1);
            upperBound = obj.fitBounds(:,2);
            options = optimoptions('fmincon','Display','notify-detailed',...
                'MaxFunctionEvaluations',10000,'MaxIterations',5000);
            fit_eq = @(b) -logLikelihood(obj.rate(b,inputVars),spike,obj.dt);
            
            [param, fval, exit_flag] = fmincon(fit_eq, initPoint,[],[],[],[],lowerBound,upperBound,[],options);
            
            if exit_flag <1 
                error(['exit flag ', num2str(exit_flag),' reached']);
            end
            
            LL = -fval;
            
        end
        

        function param = multistart(obj,data,nStarts,nTries)
            %%% Multistart fitting %%%
            % handles potential nonconvexity of optimization, minimizing
            % -LL from multiple starting points
            %
            % inputs:
            %   data - matrix of the form [timestamp, spike/no spike, input variables(1 or more)]
            %   nStarts - number of initial points for parameter estimation
            %             (parameters will estimated for each initial, this
            %             linearly increases computation time)
            %   nTries - number of times fitting is allowed to fail and restart
            %            before an error is thrown
            %
            % outputs: 
            %   param - parameters with greatest log-likelihood
            %           
            
            count = 0;
            param_options = zeros(nStarts,obj.nParams);
            LL_options = zeros(nStarts,1);
            
            for startIdx = 1:nStarts
                while true
                    try %handling if fitting doesn't work
                        initPoint = diff(obj.initBounds,1,2).*rand(obj.nParams,1)+obj.initBounds(:,1);
                        [param_options(startIdx,:),LL_options(startIdx)] = obj.fitModel(data,initPoint);
                        break
                    catch e
                        count = count+1;
                        if count == nTries
                            error(['Tried to fit this field ',num2str(nTries),' times']);
                        end
                    end
                end
            end
            
            [~,maxIdx] = max(LL_options);
            
            param = param_options(maxIdx,:);
            
        end
        
        
        function params = parameterDistribution(obj,data,iters,trainRatio)
            %%% Computing parameter distribution %%%
            % computes parameters for randomly chosen subsets of place
            % field data, dispersion of parameter values indicates
            % stability of estimates
            %
            % inputs:
            %   data - matrix of the form [timestamp, spike/no spike, input variables(1 or more)]
            %          from which subsets of data will be chosen
            %   iters - number of iterations of subsampling data and
            %           estimating parameters
            %   trainRatio - proportion of data to subsample
            %
            % outputs:
            %   params - matrix of parameter estimates, size iters x nParams
            %
            
            nDataPoints = size(data,1);
            nTrainPoints = round(trainRatio*nDataPoints);
            params = zeros(iters,obj.nParams);
            
            for itIdx = 1:iters
                
                trainInds = false(1,nDataPoints);
                shuffledInds = randperm(nDataPoints);
                trainInds(shuffledInds(1:nTrainPoints)) = true;
                trainData = data(trainInds,:);
                
                params(itIdx,:) = obj.multistart(trainData,5,5);
                
            end
            
        end
        
    end
    
end
                
                
                