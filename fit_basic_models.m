% Fit basic models to all place fields from one experimental session
%   basic models are Gaussian, Phase modulation, and PTP model
%
% inputs:
%   pathToFile - path to place field data from experimental session
%
% output:
%   no explicit outputs, parameter and log-likelihood values are saved in
%   the placefieldinfo.mat file under each place field
%

function fit_basic_models(pathToFile)

[~,filename,~] = fileparts(pathToFile);
fprintf([filename,'\n']);

load([pathToFile filesep filename '.placefieldinfo.mat'])

for fieldIdx = 1:stModel.nField
    
    fprintf('started field %1.f\n',fieldIdx)
    
    field = stModel.field{fieldIdx};
    
    % collect data
    data = vertcat(field.trial{:}); %concatenate data across trials
    temp = quantile(data(:,5),10);
    max90 = temp(end);
    data(:,5) = data(:,5)/max90; %normalize speed
    n_spikes = sum(data(:,2));
    
    % exclude fields with fewer than 100 spikes
    if n_spikes<100  
        stModel.field{fieldIdx}.gaus = [];
        stModel.field{fieldIdx}.phaseMod = [];
        stModel.field{fieldIdx}.ptp = [];
        stModel.field{fieldIdx}.baseCompare = [];
        continue
    end
    
    % make models
    gaus = gaussianModel();
    phaseMod = phaseModModel();
    ptp = ptpModel();
    modelList = {gaus,phaseMod,ptp};
    params = cell(length(modelList),1);
    
    % compute parameters for models and compare model fits
    try
        for modIdx = 1:length(modelList)
            params{modIdx} = modelList{modIdx}.parameterDistribution(data,10,.9);
        end
        
        ll = compareModels(modelList,data,50,.75);
    catch
        stModel.field{fieldIdx}.gaus = [];
        stModel.field{fieldIdx}.phaseMod = [];
        stModel.field{fieldIdx}.ptp = [];
        stModel.field{fieldIdx}.baseCompare = [];
        continue
    end
    
    stModel.field{fieldIdx}.gaus = params{1};
    stModel.field{fieldIdx}.phaseMod = params{2};
    stModel.field{fieldIdx}.ptp = params{3};
    stModel.field{fieldIdx}.baseCompare = ll;
    
end

save([pathToFile filesep filename '.placefieldinfo.mat'],'stModel')

end