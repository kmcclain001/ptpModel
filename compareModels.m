% Compare performance of place field models on data
% uses cross validation to fit each model to subsets of data and 
% test performance on remaining data  
%   
% inputs:
%   modelList - list of pfModel objects for comparison
%   data - matrix of the form [timestamp, spike/no spike, input variables(1 or more)]
%   iters - number of cross-validation iterations to perform
%   trainRatio - proportion of data to use in training each model
%                (proportion used in testing performance will be 1-trainRatio)
%
% output:
%   LLs - matrix of log-likelihoods for each model on test data, size (iters
%         x length of modelList)
%

function LLs = compareModels(modelList,data,iters,trainRatio)

nModels = length(modelList);
nDataPoints = size(data,1);
nTrainPoints = round(trainRatio*nDataPoints);
LLs = zeros(iters,nModels);


for itIdx = 1:iters
    
    trainInds = false(1,nDataPoints);
    shuffledInds = randperm(nDataPoints);
    trainInds(shuffledInds(1:nTrainPoints)) = true;
    trainData = data(trainInds,:);
    testData = data(~trainInds,:);
    
    for modIdx = 1:nModels
        
        model = modelList{modIdx};
        param = model.multistart(trainData,3);
        LLs(itIdx,modIdx) = logLikelihood(model.rate(param,testData(:,3:end)),testData(:,2),model.dt);
        
    end
    
end

end
        