%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% having run all the required sampleSizes, repeats and bootstraps on the queue 
%   (see qJob_lassoGLM(), createInputs_4_featureSelection_vs_samplesize)
% this functions evaluates the results we get from the queue_collectResults, sorting the results
% according to sampel size and repeat.
%
% INPUTS:
%   aggr_outputs: the cell returned by qOpt_collectResults_general()
%   sampleSizes : vector of sample sizes that were used to run the queue jobs
%   nRepeats    : how often did we repeat each sample size (typically 10)
%
% OUTPUTS:
%   inclusion_vs_sampleSize : nested cell array containing the inclusion probs for all repeats and sampleSize: first dimesion sample size
%   B_vs_sampleSize         : 2D-cell array storing the B-matrix returned by lassoGLM() for each sampleSize and repeat
%   fitinfo_vs_sampleSize   : 2D-cell array storing the fitinfo returned by lassoGLM() for each sampleSize and repeat

function [inclusion_vs_sampleSize, B_vs_sampleSize, fitinfo_vs_sampleSize,treeIndicesArray, successfullJobsArray] = eval_fs_vs_sampleSize_results(aggr_outputs,sampleSizes,nRepeats)

ALLOW_CONSTANT_MODEL = false; % allow a model without any variables in it?

% put the queueOutputs in more convenient format
sampleSize_tmp  = cellfun(@(a)a.sampleSize ,aggr_outputs);
repeat_tmp  = cellfun(@(a)a.repeatID ,aggr_outputs);
B_tmp = cellfun(@(a)a.B ,aggr_outputs,'UniformOutput',false);
FitInfo_tmp = cellfun(@(a)a.FitInfo ,aggr_outputs,'UniformOutput',false);
treeIndices_tmp = cellfun(@(a)a.treeIndices ,aggr_outputs,'UniformOutput',false);


% alloc our uoutput variables
inclusion_vs_sampleSize = cell(1,length(sampleSizes));
B_vs_sampleSize = cell(length(sampleSizes),nRepeats,50); %50: number of bootstraps
fitinfo_vs_sampleSize =  cell(length(sampleSizes),nRepeats,50);
treeIndicesArray =  cell(length(sampleSizes),nRepeats);
successfullJobsArray =  cell(length(sampleSizes),nRepeats);

for i = 1:length(sampleSizes) % iterate over all requested sample sizes, collect the results from the jobs that worked on this
    
    cSS = sampleSizes(i);
    ix_ss = sampleSize_tmp==cSS; %which jobs worked on that sample size
    
    inclusionProbs_across_repeats = [];
    for j = 1:nRepeats % for each repeat, get the inclusion probvability
        ix_rep = repeat_tmp == j;
        relevant_ix = find(ix_ss &ix_rep); % which jobs worked on this repeat, this sample size (typically 50, if one job per lasso call)   
        
        %for each of the N bootstraps figure out which features where included
        selectedFeatures = [];
        for k = 1:length(relevant_ix)
            tmp_ix = relevant_ix(k);
            if ALLOW_CONSTANT_MODEL
                ixLasso=B_tmp{tmp_ix}(:,FitInfo_tmp{tmp_ix}.Index1SE)~=0; % included are all that are nonZero at kappa*
            else
                ixLasso=B_tmp{tmp_ix}(:,FitInfo_tmp{tmp_ix}.Index1SE)~=0; % included are all that are nonZero at kappa*
                if all(ixLasso==0) % an all constant model! use the next complex model
                    assert(FitInfo_tmp{tmp_ix}.Index1SE == size(B_tmp{tmp_ix},2),'if all features are not important, the index shuold be the last')
                    ixLasso = B_tmp{tmp_ix}(:,end-1)~=0;
                    warning('constant model not allowed')
                end
            end
            selectedFeatures=[selectedFeatures; ixLasso'];
            
            %additionally, remember the lasso fit for later
            B_vs_sampleSize{i,j,k} = B_tmp{tmp_ix};
            fitinfo_vs_sampleSize{i,j,k} = FitInfo_tmp{tmp_ix};
            
        end
        tmp_inclusion = mean(selectedFeatures); %the inclusion prob is just the average over the 0/1 across bootstraps
        treeIndicesArray{i,j} = treeIndices_tmp(relevant_ix);
        successfullJobsArray{i,j} = length(relevant_ix);
        
        %no single job returned stuff
        if isempty(relevant_ix)
            tmp_inclusion = nan(1,size(inclusionProbs_across_repeats,2));
        end
        
        inclusionProbs_across_repeats = [inclusionProbs_across_repeats; tmp_inclusion];
    end
    inclusion_vs_sampleSize{i} = inclusionProbs_across_repeats;
end
