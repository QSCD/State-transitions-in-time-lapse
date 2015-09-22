%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
%% Figure 5A,B,C,D
% given the results from the power-simulations (stored in powerTophat)
% create the plots used in Figure 5 of the paper

load data/power/powerTophat.mat
addpath('..')
trueFeaturesFlag = [0 0 0 0  ones(1,7)  0 0 0 ];

% power
powerPlot(powerArray,sampleSizes,effectStrength,featureNames )


% enlarged plot for phi_2
phi2_ix = 7;
powerPlot(powerArray(:,:,phi2_ix),sampleSizes,effectStrength,featureNames(phi2_ix) )
hline(0,'r') %indicate where the sample size plot in Fig4 comes from
colorbar


cSS = 4500;   % same setting as for Fig4
cEffectSize = 1; % same setting as for Fig4
ii = effectStrength == cEffectSize;
sampleSizePlot(sampleSizes, squeeze(powerArray(ii,:,:)), trueFeaturesFlag)
set(gca,'XTick',[0:1000:4000])
set(gca,'YTick',[0:0.2:1])
ylabel('Power')

%separate the relevant and irreleavnt features for the powerPlots
%relevant features
subplotLayout = [3,3];
powerPlot(powerArray(:,:,trueFeaturesFlag==1),sampleSizes,effectStrength,featureNames(trueFeaturesFlag==1), subplotLayout )

%irrelevant features
subplotLayout = [3,3];
powerPlot(powerArray(:,:,trueFeaturesFlag==0),sampleSizes,effectStrength,featureNames(trueFeaturesFlag==0), subplotLayout)
