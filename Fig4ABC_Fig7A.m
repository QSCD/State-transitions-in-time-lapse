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

% transform trees into feature matrix
if exist('data/tophat/dataset_tophatDependence.mat','file')~=2
    createFeatureMatrix('data/tophat/trees_tophatDependence.mat', 'data/tophat/dataset_tophatDependence','marker')
end

%% Fig4
%==================
%load full data
load data/tophat/dataset_tophatDependence.mat

fulldata = [X;X_not];
fulllabels = [ones(length(X),1); zeros(length(X_not),1)];

%==================
% do lasso several times
% this can take quite a while, so take the precalculated results if available
if exist('data/tophat/featureSelection_tophatDependence.mat','file')~=2
    
    [selectedFeatures, B, fitinfo] = lassoFeatureSelectionWithReplication(zscore(fulldata ),fulllabels,attributes,50,'poisson');
    save data/tophat/featureSelection_tophatDependence.mat selectedFeatures B fitinfo fulldata attributes
end

%==================
% plotting the bootstrapped results
trueFeaturesFlag = [0 0 0 0  ones(1,7)  0 0 0 ];
plotResults_LassoBootstrap_matfile('data/tophat/featureSelection_tophatDependence.mat',0,trueFeaturesFlag)


%=======================
% debiasing + extracting the true weights of the hazard
load data/tophat/featureSelection_tophatDependence.mat
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
%plot the true kernel contribution
tmpx = -150:150;
trueKernel = (abs(tmpx)<70).*0.0025;
hold on
plot(tmpx,trueKernel,'k-')
xlim([-100 100])
ylabel('Density weights')
xlabel('Distance r')
box off



%% Fig7 : sister correlations 
load data/tophat/trees_tophatDependence.mat trees_notPruned_withDensities
%NOTICE: THESE FUNCTION WAS DEISGN TO HANDLE THE WEIGHTS THE COME OUT OF LASSO, hence the expect the
%minus sing and mulitplied with dt!!
b_helper = -b_true.*dt;
[pval_true,O_true,E_true] = checkIndependenceOfSisters_fullModel(trees_notPruned_withDensities, 'marker', b_helper, selectedAttributes, NOCONSTANT_FLAG);


% -----------------------------------
% deliberately remove density functions: the IMPAIRED scenario
load data/tophat/trees_tophatDependence.mat trees_notPruned_withDensities
load data/tophat/featureSelection_tophatDependence.mat

wrong_selected = selectedFeatures;
wrong_selected(:,10:11)=0;

dt = 0.1;
NOCONSTANT_FLAG = 1;
[wrong_b,wrong_selectedAttributes] = extractWeigh_LassoBootstrap(wrong_selected,attributes,fulldata,fulllabels,dt,NOCONSTANT_FLAG);
wrong_b_helper = -wrong_b.*dt;
[pval_impaired,O_impaired,E_impaired] = checkIndependenceOfSisters_fullModel(trees_notPruned_withDensities,'marker',wrong_b_helper,wrong_selectedAttributes,NOCONSTANT_FLAG);

% ---------------------------------
% put the above into the same plot!
assert(all(O_true == O_impaired),'observations should be the same');

plot_sister_correlation_figure(O_true,E_true,E_impaired, pval_true,pval_impaired  )
printfig('sisterCorrStuff','pdf', 'figfile',true,  'FontSize',4, 'LineWidth', 0.5, 'FigureWidth', 'halfpage')

