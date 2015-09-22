%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

%plots the estimated differentiationRate including confidence intervals vs the true one
% indepVar : the variable for which we calculate the diff.rate (eg, time or density)
% MAP      : the estimated hazard (events/all)
% lowerCI,upperCI: confidence region (in case of count data, look at estimateProportions.m, it uses a beta distribtution to get confidence intervals)
function differentiationRatePlot(indepVar,MAP,lowerCI,upperCI,truth)

    assert(length(indepVar)==length(MAP));
    assert(length(indepVar)==length(lowerCI));
    assert(length(indepVar)==length(upperCI));
    
    if ~isempty(truth)
        assert(length(indepVar)==length(truth));
    end

    upperError = upperCI-MAP; % <- because shadedErrorbar expects the magnitude of the error not the abs position
    lowerError = MAP-lowerCI; % <- because shadedErrorbar expects the magnitude of the error not the abs position
    handles = shadedErrorBar(indepVar,MAP,[upperError;lowerError],[],false,false);
    hold on

    if ~isempty(truth)
        handles.truth = plot(indepVar,truth,'-r');
        legend([handles.mainLine,handles.patch, handles.truth],{'Posterior mean \hat\lambda','95% CI','true'},'Location','NorthWest')
    else
        legend([handles.mainLine,handles.patch],{'Posterior mean \hat\lambda','95% CI'},'Location','NorthWest')
    end
    ylabel('Differentiation rate (1/h)');
    
    box off
end
