%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

%% creates Fig3 CDE of the paper
addpath(genpath('algorithm'))
addpath(genpath('utils'))
addpath(genpath('plotting'))

dt = 0.1; % the observation interval is 0.1. important for scalign the lambda correctly: lambda*dt = #events/#total

% loading + pruning + removing duplicates
load data/tophat/trees_tophatDependence.mat trees_notPruned_withDensities
pC_time = @(t)t.('marker')==1;g = @(mT)tUtil_pruneMichiTree_TIME(mT,pC_time); %prune the trees at the differentiation event
treesPrunedwithDensity = cellfun(g, trees_notPruned_withDensities,'UniformOutput',0);clear trees_notPruned_withDensities;
treesPrunedwithDensity = cellfun(@(t)removeDuplicateTimepoints(t),treesPrunedwithDensity,'UniformOutput',0);


%% Figure 3C: Density dependence dendence of lambda
% ------ tophat 70 --------
% kernel with radius 70
nBins = 47;
tophat70_did = cellfun(@(t)t.density_70(t.marker==1),treesPrunedwithDensity,'UniformOutput',false);    tophat70_did = [tophat70_did{:}];
tophat70_didnot = cellfun(@(t)t.density_70(t.marker==0),treesPrunedwithDensity,'UniformOutput',false);tophat70_didnot = [tophat70_didnot{:}];


[~, ~, dnes70range, lambda_dens70, lowerCI_dens70, upperCI_dens70] = estimateDifferentiationRate(tophat70_did,tophat70_didnot,nBins,dt);
trueDens70_hazard = 0.0025 .* dnes70range;
figure;
differentiationRatePlot(dnes70range,lambda_dens70,lowerCI_dens70,upperCI_dens70,trueDens70_hazard);xlabel('Density \rho_{70}')
xlim([0 40])
ylim([0 0.12])



%% Figure 3C inset: Time dendence of lambda
nBins=50;

time_did = cellfun(@(t)t.absoluteTime(t.marker==1),treesPrunedwithDensity,'UniformOutput',false);    time_did = [time_did{:}];
time_didnot = cellfun(@(t)t.absoluteTime(t.marker==0),treesPrunedwithDensity,'UniformOutput',false);time_didnot = [time_didnot{:}];

[~, ~, timerange, lambda_time, lowerCI_time, upperCI_time] = estimateDifferentiationRate(time_did,time_didnot,nBins,dt);
figure;
differentiationRatePlot(timerange,lambda_time,lowerCI_time,upperCI_time,[]);xlabel('Time (h)')
xlim([0,100])
ylim([0,0.035])


%% Figure 3D: Time dendence + density dependence of lambda

t_tmp = sort(unique([time_did time_didnot]));
d_tmp = sort(unique( [tophat70_did tophat70_didnot]));
edges = {linspace(min(t_tmp),max(t_tmp),30)   ,   linspace(min(d_tmp),max(d_tmp),30)};

n =hist3([time_did' tophat70_did'],'Edges',edges);
n_Back = hist3([[time_did time_didnot]'  [tophat70_did tophat70_didnot]'],'Edges',edges);

lambda_time_dens = (n./n_Back)./dt;
lambda_time_dens(isnan(lambda_time_dens))=0;


figure;
surf(edges{1},edges{2},log10(lambda_time_dens'),'LineStyle','none')
ylabel('Density \rho_{70}');
xlabel('Time')
xlim([0 97])
ylim([0 38])
view([0 90])
grid off
colorbar
