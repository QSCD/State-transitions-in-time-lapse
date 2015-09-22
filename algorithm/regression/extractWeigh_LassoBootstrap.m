%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% having done the feature selection, we now refit the model with the selected features
% -> DEBIASING
%
% also we transform back the estimated weights into the true weights of the hazard wTrue = -w./dt
%
% INPUT
%   selectedFeatures: binary matrix of selected features across bootstraps (n_bootstraps x features)
%   attributes      : cell array of attribute names  (features x 1
%   fulldata        : full datamatrix (samples x features)
%   fulllabels      : full class label vector (samples x 1)
%   dt              : observation interval, used to rescale the weights
%   
% OUTPUT
%   newweights      : fitted, unbiased weights vector of the selected attributes
%   newattributes   : names of the selected attributes
function [newweights,newattributes] = extractWeigh_LassoBootstrap(selectedFeatures,attributes,fulldata,fulllabels,dt,NOCONSTANT_FLAG)

INLCUSION_THRESHOLD = 0.9;

%which features are about the threshold -> include them in the final model
ixLasso_final  = mean(selectedFeatures,1)>=INLCUSION_THRESHOLD;
disp('included features:')
attributes(ixLasso_final)

% Lasso was fitted with an offset term (which cannot be penalized)
% should the debiased regression be with or without the offset
% actually, doesnt make a difference
if NOCONSTANT_FLAG == 0
    [b_all] = glmfit(fulldata(:,ixLasso_final),~fulllabels, 'poisson');
else
    [b_all,dev,stats]= glmfit(fulldata(:,ixLasso_final),~fulllabels, 'poisson','constant','off');
end

%since we estimated not the weights directly but instead w = -trueW.*dt *X
% we get back trueW
newweights = -b_all./dt;
newattributes = attributes(ixLasso_final);
