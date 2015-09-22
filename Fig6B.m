%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% Figure 6B
% How does tracking error impact on the selection of the features

nRepeats = 10;
trackingErrorVector = sort([0.001 0.005 0.01 0.015 0.02 0.05 0 0.1]);

load data/trackingError/trackingError_gaussTime;

% count how often the feature ewas included, for each tracking error amount
[inclusion_vs_tracking,Barray,fitArray,tIndexArray] = eval_fs_vs_sampleSize_results(aggr_outputs,trackingErrorVector,nRepeats);

trueFeaturesFlag = [1 0 0 0 ones(1,5) 0 0  0 0 0]; % time and a gaussian kernel

%-----------------------
%plot the power
power = cellfun(@(t)mean(t>=0.9),inclusion_vs_tracking,'UniformOutput',false); power = vertcat(power{:});

%slightly abusing this sample size plotting function, works also for
%tracking error on the x-axis

sampleSizePlot(100.*trackingErrorVector, power, trueFeaturesFlag)
xlabel('Tracking error (fraction of cells) in percent');
ylabel('Power');
set(gca,'XTick',[0:1:10])
set(gca,'YTick',[0:0.2:1])

