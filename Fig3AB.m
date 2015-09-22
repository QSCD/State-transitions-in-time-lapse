%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

%% Figure 3A,B
addpath(genpath('algorithm'))
addpath(genpath('utils'))
addpath(genpath('plotting'))
%--------------------
% constant time
clear
dt = 0.1;
c=0.01;  %true constant rate

%load data
load data/constant_time/trees_constant.mat
trees_const = trees_notPruned_withDensities;
pC_time = @(t)t.('marker')==1;g = @(mT)tUtil_pruneMichiTree_TIME(mT,pC_time);
trees_const = cellfun(g, trees_const,'UniformOutput',0);
trees_const = cellfun(@(t)removeDuplicateTimepoints(t),trees_const,'UniformOutput',0);

%do the estimate of lambda
nBins=20;
time_did = cellfun(@(t)t.absoluteTime(t.marker==1),trees_const,'UniformOutput',false);
time_did = [time_did{:}];
time_didnot = cellfun(@(t)t.absoluteTime(t.marker==0),trees_const,'UniformOutput',false);
time_didnot = [time_didnot{:}];
[~, ~, timerange, lambda_time, lowerCI_time, upperCI_time] =estimateDifferentiationRate(time_did,time_didnot,nBins,dt);

%plotting
figure;differentiationRatePlot(timerange,lambda_time,lowerCI_time,upperCI_time,c.*ones(size(timerange)));xlabel('Time (h)')
xlim([0,100])
ylim([0,0.03])



%--------------------
% linear time
clear

dt = 0.1;
slope =3e-4;  % true slope of the differentiation rate

%load data
load data/constant_time/trees_LinearTime
trees_linear = trees_notPruned_withDensities;
pC_time = @(t)t.('marker')==1;g = @(mT)tUtil_pruneMichiTree_TIME(mT,pC_time);
trees_linear = cellfun(g, trees_linear,'UniformOutput',0);
trees_linear = cellfun(@(t)removeDuplicateTimepoints(t),trees_linear,'UniformOutput',0);

%do the estimate of lambda
nBins=20;
time_did = cellfun(@(t)t.absoluteTime(t.marker==1),trees_linear,'UniformOutput',false);
time_did = [time_did{:}];
time_didnot = cellfun(@(t)t.absoluteTime(t.marker==0),trees_linear,'UniformOutput',false);
time_didnot = [time_didnot{:}];
[~, ~, timerange, lambda_time, lowerCI_time, upperCI_time] =estimateDifferentiationRate(time_did,time_didnot,nBins,dt);

%plotting
figure;differentiationRatePlot(timerange,lambda_time,lowerCI_time,upperCI_time,slope.*timerange);xlabel('Time (h)')
xlim([0,100])
ylim([0,0.03])
