%% Demo for Position Theta Phase model 
%
% Author: Kathryn McClain
% Email: km3911@nyu.edu
% Date: 7/23/19
%
% Model presented in McClain et al. 2019

%% load data
load('example_dataset_DT2_20160227.placefieldinfo.mat');

%choose example place field
field = stModel.field{150} %stModel is struct with place field data

% combine data across trials
data = vertcat(field.trial{:});

%% plot place field summary info
figure();hold on
subplot(2,1,1)
fr = tuning_curve(data(:,2),data(:,3),linspace(0,1,51),1/1250);
plot(linspace(0,1,50),fr);
ylabel('Trial-averaged firing rate')
xlabel('Position in field')
fix_ax(gca)

subplot(2,1,2)
scatter(data(data(:,2)>0,3),data(data(:,2)>0,4),'.');
xlabel('Position in field')
ylabel('Theta phase')
fix_ax(gca)

%% Fit PTP model to palce field data

% create PTP model object
ptp = ptpModel();

% compute parameters for 5 random subsets of data
param_set = ptp.parameterDistribution(data,5,.9);

% estimate parameter as median parameter value
b = median(param_set);

%% Plot PTP model

figure();hold on

subplot(2,2,1)
x = linspace(0,1,1000);
phs = linspace(0,2*pi,1000);
[X,Phs] = meshgrid(x,phs);
R = exp(b(1)+((-1/2)*((X-b(2))/b(3)).^2)).*...
    exp(b(4)*(cos(Phs-(b(6)+ b(5)*(X-b(2))))-1));
im = imagesc(x,phs,R);
colormap('summer')
set(gca,'YDir','normal');
im.AlphaData = 1;
ylim([0,2*pi]);
fix_ax(gca)
ylabel('Theta phase')
title('PTP model')
h = colorbar;
ylabel(h,'Firing rate')
fix_ax(h)

subplot(2,2,2)
thetaModFunc = exp(b(4)*(cos(phs-b(6))-1)); % at middle of place field
plot(phs,thetaModFunc);
xlabel('Theta phase (rad)')
ylabel('Phase preference')
title('Theta Modulation') 
fix_ax(gca)

subplot(2,2,3)
spatialFunc = exp(b(1)+((-1/2)*((x-b(2))/b(3)).^2));
plot(x,spatialFunc);
xlabel('Position in field');
ylabel('Firing rate')
title('Spatial input')
fix_ax(gca)


subplot(2,2,4)
precessionFunc = b(6)+ b(5)*(x-b(2));
plot(x,precessionFunc);
xlabel('Position in field')
ylabel('Preferred phase')
title('Precession')
fix_ax(gca)

%% Simulate experiments with model fit to example place field
% two options for simulation:
%   1) simulate from scratch, i.e. simulate position and theta phase then
%      simulate firing rate and spiking from simulated inputs
%   2) use real measured input variables and simulate firing rate and spiking
%      from these
%

% 1) simulate input variables for experiment
inputVars = simulateTrial(.6,8.8,0,30);

% compute firing rate with PTP model and parameters from example place
% field
firing_rate = ptp.rate(b,inputVars);

% simulate spikes through poisson process
spikes = poissrnd(firing_rate.*ptp.dt); 

% 2) alternatively, simulate using input variables measured in real
% experiment
firing_rate_alt = ptp.rate(b,data(:,3:4));
spikes_alt = poissrnd(firing_rate_alt.*ptp.dt);

%% Plot simulated data

figure();hold on
subplot(2,2,1);hold on
cross = [0;find(diff(inputVars(:,2))<0);size(inputVars,1)];
for i = 1:(length(cross)-1)
    plot_inds = (cross(i)+1):cross(i+1);
    plot(inputVars(plot_inds,1),inputVars(plot_inds,2),'k');
end
xlabel('Position in field')
ylabel('Theta phase (rad)')
title('Simulated trial')
fix_ax(gca)

subplot(2,2,2); hold on
plot(inputVars(:,1),firing_rate,'b');
scatter(inputVars(spikes>0,1),2*ones(sum(spikes>0),1),'k.');
xlabel('Position in field')
ylabel('Firing rate (Hz)')
title('Simulated trial')
fix_ax(gca)

subplot(2,2,3)
fr = tuning_curve(spikes_alt,data(:,3),linspace(0,1,51),1/1250);
plot(linspace(0,1,50),fr);
ylabel('Trial-averaged firing rate')
xlabel('Position in field')
title('Simulated experiment')
fix_ax(gca)

subplot(2,2,4)
scatter(data(spikes_alt>0,3),data(spikes_alt>0,4),'.');
xlabel('Position in field')
ylabel('Theta phase')
title('Simulated experiment')
fix_ax(gca)

%% Compare basic models (as in Figure S3)

% create model variants
gaus = gaussianModel();
phaseMod = phaseModModel();

% compute parameters for variants
param_set_gaus = gaus.parameterDistribution(data,5,.9);
param_set_phaseMod = phaseMod.parameterDistribution(data,5,.9);
b_gaus = median(param_set_gaus);
b_phaseMod = median(param_set_phaseMod);

% compute cross-validated LLs
modelList = {gaus,phaseMod,ptp};
LL = compareModels(modelList,data,10,.75);

%% Plot model comparison

figure();hold on

subplot(2,1,1); hold on
gaus_fr = gaus.rate(b_gaus,inputVars);
phaseMod_fr = phaseMod.rate(b_phaseMod,inputVars);
plot(inputVars(:,1),gaus_fr,'r');
plot(inputVars(:,1),phaseMod_fr,'c');
plot(inputVars(:,1),firing_rate,'b');
legend('Gaussian model','Phase mod','PTP')
xlabel('Position in field')
ylabel('Firing rate')
title('Firing rate predicted for 3 models in simulated trial')
fix_ax(gca)

subplot(2,1,2);
c = categorical({'Gaussian','Phase mod','PTP'});
norm_ll = zscore(mean(LL));
bar(c,norm_ll);
ylabel('Zscore log-likelihood')
title('Comparison of 3 models performance in real data')
fix_ax(gca)
