%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% does GLM + Lasso to obtain the relevant features
% resamples the dataset "N_replications"-times
%
% INPUT
%   fullData        : datamatrix:  samples x features
%   fullLabels      : vector (samples x 1) of class labels (0= undiff, 1=differentiating)
%   featureNames    : cell array (features x 1) of feature names 
%   N_replications  : number of "bootstraps" to perform 
%
% OUTPUT
%   selectedFeatures        : binary matrix of selected features in each bootstrap run (N_replications x features)
%   B_array,FitInfo_array   : cell array (N_replications x 1) containing the info provided by lassoGLM
%   
function [selectedFeatures,B_array,FitInfo_array] =  lassoFeatureSelectionWithReplication(fullData,fullLabels,featureNames,N_replications,regressionType)

assert(length(featureNames)==size(fullData,2));

selectedFeatures = zeros(N_replications,length(featureNames)); %for each "bootstrap",its 1 if feature was selected or 0 if not
B_array = cell(1,N_replications);
FitInfo_array = cell(1,N_replications);

FULLBOOTSTRAP_FLAG = 1;
parfor i = 1:N_replications
    disp(i)
    
    %subsample the majority class
    ix_nonEvent = find(fullLabels==0);
    ix_Event = find(fullLabels==1);
    
    tmp_ix = unidrnd(length(ix_nonEvent),3.*length(ix_Event),1);
    subsample_ix_noEvent = ix_nonEvent(tmp_ix);
    assert(all(fullLabels(subsample_ix_noEvent)==0))
    
    %bootstrap/resampel also the minority class!
    if FULLBOOTSTRAP_FLAG
        tmp_ix_events = unidrnd(length(ix_Event),length(ix_Event),1);
        subsample_ix_events = ix_Event(tmp_ix_events); 
    %otherwise just take all events are they are
    else 
        subsample_ix_events = ix_Event;
    end
    
    subsample_data = [fullData(subsample_ix_events,:); fullData(subsample_ix_noEvent,:)];
    subsample_label= [fullLabels(subsample_ix_events); fullLabels(subsample_ix_noEvent)];

    %% do the feature selection
    
    %for poisson reg, invert the labels
    if strcmp(regressionType,'poisson')
        subsample_label = double(~subsample_label);
    end
    [B,FitInfo] = lassoglm(subsample_data,subsample_label,regressionType,'CV',10,'PredictorNames',featureNames);
    
    % take either the features at crossval lambda (if that model is not the constant)
    % or (if the cross val model is constant), the first non-constant model!
    
    if all(B(:,FitInfo.Index1SE) == 0)
        FitInfo.Index1SE = FitInfo.Index1SE-1; % take the model which is 
        FitInfo.nonConstantModelEnforced = true;
    else
        FitInfo.nonConstantModelEnforced = false;
    end
    ixLasso=B(:,FitInfo.Index1SE) ~= 0;
    selectedFeatures(i,:)=ixLasso;
    B_array{i} = B;
    FitInfo_array{i} = FitInfo;
    
end
