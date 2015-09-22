%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

addpath(genpath('algorithm'))
addpath(genpath('utils'))
addpath(genpath('plotting'))

if exist('data/gaussian_time/dataset_Gaussian_timeDependence.mat','file')~=2
    createFeatureMatrix('data/gaussian_time/trees_Gaussian_timeDependence.mat', 'data/gaussian_time/dataset_Gaussian_timeDependence', 'marker')
end

%% GLM
% load and create some data
load data/gaussian_time/dataset_Gaussian_timeDependence.mat

fulldata = [X;X_not];
fulllabels = [ones(length(X),1); zeros(length(X_not),1)];

%==================
% do lasso several times
% this can take quite a while, so take the precalculated results if available
if exist('data/gaussian_time/featureSelection_Gaussian_timeDependence.mat','file')~=2
    [selectedFeatures, B, fitinfo] = lassoFeatureSelectionWithReplication(zscore(fulldata),fulllabels,attributes,50,'poisson');
    save data/gaussian_time/featureSelection_Gaussian_timeDependence.mat selectedFeatures B fitinfo attributes
end

%==================
% plotting
trueFeaturesFlag = [1 0 0 0 ones(1,5) 0 0 0 0 0];
plotResults_LassoBootstrap_matfile('data/gaussian_time/featureSelection_Gaussian_timeDependence.mat',0,trueFeaturesFlag)


%=======================
% debiasing + extracting the true weights of the hazard
load data/gaussian_time/featureSelection_Gaussian_timeDependence.mat
dt = 0.1;
NOCONSTANT_FLAG = 1;
[b_true,selectedAttributes] = extractWeigh_LassoBootstrap(selectedFeatures,attributes,fulldata,fulllabels,dt,NOCONSTANT_FLAG);

%==================
% the kernel estimated from the regression
% --------------------
% having done the above feature selection with the lasso and replication,
%  fit the sparse model to the full data one last time (debiasing)

kernelWeights = extractDensityFeatures(b_true,attributes,selectedAttributes,NOCONSTANT_FLAG);
reconstruct_Densitykernel_from_rings(kernelWeights)

% --------------------
%plot the true kernel contribution
tmpx = -150:150;

trueKernel = normpdf(tmpx,0,30);
trueKernel = trueKernel./sum(trueKernel);

hold on
plot(tmpx,trueKernel.*0.38,'k-') % since the kernel was not normalized in the simulation
xlim([-100 100])
box off
ylabel('Density weights')
xlabel('Distance r')
