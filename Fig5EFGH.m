%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% given the results from the power-simulations (stored in powerGaussTime_xxx.mat)
% create the plots used in Figure 5 of the paper
%
% note that parameters for time and density are varied seperately, hence we
% have to blend those two
addpath('..')

q = load('data/power/powerGaussTime_densityVariation'); % created in 'newPower_varyingDensity'
w = load('data/power/powerGaussTime_TimeVariation'); % created in 'newPower_varyingTime'

trueFeaturesFlag = [1 0 0 0  ones(1,5) 0 0 0 0 0 ];

%plot them seperately
powerPlot(q.powerArray,q.sampleSizes,q.effectStrength,q.featureNames )
powerPlot(w.powerArray,w.sampleSizes,w.effectStrength,w.featureNames )


% when looking at the powerplot of the time variable, it makes
% only sense if the time parameter is varied, similarly for the density
% hence merge those plots, combining those where the features are actually
% the ones that are varied

mergedPower = q.powerArray;
mergedPower(:,:,1) = w.powerArray(:,:,1);
mergedPower(:,:,2) = w.powerArray(:,:,2);

%relevant features
subplotLayout = [3,3];
powerPlot(mergedPower(:,:,trueFeaturesFlag==1),q.sampleSizes,q.effectStrength,q.featureNames(trueFeaturesFlag==1), subplotLayout )

%irrelevant features
subplotLayout = [3,3];
powerPlot(mergedPower(:,:,trueFeaturesFlag==0),q.sampleSizes,q.effectStrength,q.featureNames(trueFeaturesFlag==0), subplotLayout)


cSS = 4500;   % same setting as for Fig4
cEffectSize = 1; % same setting as for Fig4
ii = q.effectStrength == cEffectSize;
sampleSizePlot(q.sampleSizes, squeeze(mergedPower(ii,:,:)), trueFeaturesFlag)
set(gca,'XTick',[0:1000:4000])
set(gca,'YTick',[0:0.2:1])
ylabel('Power')


