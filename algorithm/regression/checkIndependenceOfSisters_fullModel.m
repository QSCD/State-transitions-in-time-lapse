%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

%
% given some trees and a function that gives us probabilities for cells to diff (dep. on exteral factor)
% how many differentiation events of sister pairs do we expect and how many do we observe
% are they compatible?
%
%   trees_with_density: just some trees with all neccesary features
%   marker : field of the tree that tells us wether the cell is undiff/diff
%   beta   : weights from the GLM
%   correspondingAttributes: the features corresponding to the betas
%   NOconstantTermFlag: allow an offet w0?
%
% returns the pvalue of the chi2 test, and the observed, expected frequencies of the trees
function [pval, Observed,Expected] = checkIndependenceOfSisters_fullModel(trees_with_density,marker,beta,correspondingAttributes,NOconstantTermFlag)

sisterProbs = []; % [p1 p2] for each sister pair in the data: prob to differntiate for each sister cell given the external hazard
eventHappend = []; %  did both sisters diff or one or none

%if not done already put the rings as features
if ~isfield(trees_with_density{1},'ring_density_10')
    trees_with_density = add_RingDensityFeatures(trees_with_density);
end

for i=1:length(trees_with_density)
   
    cT = trees_with_density{i}; 
    cells = unique(cT.cellNr);  

    %we only have to look at cells that are undifferentiated completely
    % only those can have two daughter with onsets (or one onset, or no onset)
    relevantCells= findCompleteUndiffCells(cT,marker);
    for j=1:length(relevantCells)
        cC = relevantCells(j);
        
        %check if the cell has two daughters and is not differentiated itself
        lD = cC.*2;
        rD=lD+1;
        
        if all(cT.(marker)(cT.cellNr==cC)==0) % the cell itself is undifferentiated
           if any(cells == lD) & any(cells == rD) % and the two daughter exist
                             
                %was there an event in any of the two cells?
                e1 =any(cT.(marker)(cT.cellNr==lD));
                e2 =any(cT.(marker)(cT.cellNr==rD));
                eventHappend = [eventHappend; e1 e2];
               
                % make the two cells into a data matrix
                datMat1 = cell2dataMatrix(cT,lD,correspondingAttributes);
                datMat2 = cell2dataMatrix(cT,rD,correspondingAttributes);
                p1 = calc_prob2diff(beta,datMat1,NOconstantTermFlag);
                p2 = calc_prob2diff(beta,datMat2,NOconstantTermFlag);
                sisterProbs = [sisterProbs; p1 p2];

           end
        else 
            error('should never get here')
            %if the cC is nowhere undiffereniated, theres no need to consider daughters, so mark themn to be skipped
            subtree = tUtil_getSubtree(cT,cC) ;
           cells2skip = [cells2skip unique(subtree.cellNr)];
           
        end
    end
end

%% do the chi2test to check if observed and expected proportions are equal
% are the observed pairs statistically compatible with what we expect?
[pval, Observed, Expected] = doIndependenceAnalysis(eventHappend,sisterProbs(:,1),sisterProbs(:,2));

%plot the bars seperately
figure;
colormap([0 0 0;0.7 0.7 0.7])
bar([Observed;Expected]')
set(gca,'XTickLabel',{'both','one','none'})
legend({'Observed','Expected'})
box off
ylabel('absolute frequency')
title(num2str(pval))
end


% calcualtes the probability that the cell with cellID differentiated, given some estimated  rate
% we plug in the current "environment" of the cell into the lamda for that
function p= calc_prob2diff(beta,datMat,NOconstantTermFlag)
    beta = -beta;
    if NOconstantTermFlag == 0
        if length(beta)==1 %its a model with only the constant term!
            theRate = beta(1);
        else               %its a model with constant term + some others
            theRate = beta(1) + datMat*beta(2:end);
        end
    else                   %its a model withOUT constant term
        theRate = datMat*beta(1:end);
    end
    if any(isnan(theRate)) | any(theRate<0)
        warning('some hazard is nan, possibly because extrapolation')
    end
    theRate(isnan(theRate))=0;

    % dt is already absorbed into the weights!!
    % no need to consider it explicitly
    % therefore:
    % prob to diff is 1-exp(-int rate)
    % prob to not diff is exp(-int rate) 

    % we want the prob to not diff in each observed interval and multiply them over the whole cell =
    % exp(-i1) * exp(-i2 ) ... = exp(-sum(i))
    % the prob to differentiate anywhere inside this cell is 1-P(not diff inside the cell) = 1-exp(-sum(i))
    p = 1-exp(-sum(theRate));
    assert(~isnan(p),'returns NaN!!')
end

%finds all cells that are completely undifferentiated
function relevantCells= findCompleteUndiffCells(tree,marker)
    allcells = unique(tree.cellNr);
    relevantCells = [];
    for i = 1:length(allcells)
        cC = allcells(i);
        if all(tree.(marker)(tree.cellNr==cC)==0)
            relevantCells = [relevantCells cC];
        end
    end
end


%% given a tree, a cellNr and some features, make a dataMatrix (rows=samples, cols = features) out of it
function datMat = cell2dataMatrix(cT,cCell,correspondingAttributes)

    ix = cT.cellNr ==cCell;
    datMat = zeros(sum(ix),length(correspondingAttributes));

    for i  =1:length(correspondingAttributes)
        datMat(:,i) = cT.(correspondingAttributes{i})(ix);
    end
end

% calculates the pval of the observed pairings under independence assumption
%
%  pairs: Nx2 array of zeros and ones: wether one both or no sister differentiated
%  p1: Nx1 array of marginal prob that sister on differentiates (determined via its environment eg)
%  p2: Nx1 array: analogous to p1
function [pval, Observed, Expected] = doIndependenceAnalysis(pairs,p1,p2)

    % what have we observed
    tmp_o = sum(pairs,2);
    both_o = sum(tmp_o==2);
    one_o = sum(tmp_o==1);
    none_o = sum(tmp_o==0);
    Observed = [both_o one_o none_o];

    % what do we expect in terms of two one none pairs
    both_e = sum(p1.*p2);
    one_e = sum(p1.*(1-p2) + p2.*(1-p1));
    none_e = sum((1-p1).*(1-p2));
    Expected = [both_e one_e none_e];

    % do the chi^2 test (two degrees of freedom)
    C = sum( ((Observed-Expected).^2)./Expected );
    pval = 1-chi2cdf(C,2);

end