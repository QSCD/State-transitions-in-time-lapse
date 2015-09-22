%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% INPUT:
%   events: the atributes of the event (e.g at what time or density it happend)
%   nonEvents: the same attrubte for all datapoint where the event did not happen
%   dT: the observatino period of the sample
%
% OUTPUT:
%   n,n_background, range: the hisogram data
%   MAP                  : posterior mea estimate
%   lowerCI              : lower boundary of cred. interval
%   upperCI              : upper boundary of cred. interval
function [n, n_background, range, MAP, lowerCI, upperCI] = estimateDifferentiationRate(events,nonEvents,nBins,dt)

    if ~iscell(events) & ~iscell(nonEvents) % if its no cell array, create a one element array
        events = {events};
        nonEvents = {nonEvents};
    elseif iscell(events) & iscell(nonEvents)
        assert(length(events)== length(nonEvents),'events/nonevents has to have same size')
    else
       error('both events and nonEvents must be arrays or both must be non arrays') 
    end

    for i = 1:length(events)

        unionOfBoth= [events{i} nonEvents{i}];
        range = linspace(min(unionOfBoth),max(unionOfBoth),nBins);
        [n_background, ~] = hist(unionOfBoth,range); % histogram of background 
        [n, ~] = hist(events{i},range); % histogram of the events
        %get some confidence intervals
        [MAP, lowerCI, upperCI]= estimateProportions(n,n_background);
    end

    MAP=MAP./dt;
    lowerCI=lowerCI./dt;
    upperCI=upperCI./dt;
end

%estimate the fraction and its credibility intervals
% Jeffreys Prior -> alpha = beta = 1/2
% uniform prior  -> alpha = beta = 1
function [MAP, lowerCI, upperCI]= estimateProportions(nSucesses,nOverall)
    assert(length(nSucesses) == length(nOverall))
    prior_constant = 0.5;

    CI_alphaLevel = 0.05; % 95% CI
    lowerCI = zeros(1,length(nSucesses));
    upperCI = zeros(1,length(nSucesses));

    alphasArray= zeros(1,length(nSucesses));
    betasArray= zeros(1,length(nSucesses));
    for i = 1:length(nSucesses)
        posteriorAlpha= prior_constant+nSucesses(i);
        posteriorBeta= prior_constant+nOverall(i)-nSucesses(i);

        q95 = betainv(1-CI_alphaLevel,posteriorAlpha,posteriorBeta);
        q5 = betainv(CI_alphaLevel,posteriorAlpha,posteriorBeta);

        lowerCI(i) =q5;
        upperCI(i) =q95;

        alphasArray(i) =posteriorAlpha;
        betasArray(i) =posteriorBeta;
    end

    MAP = (nSucesses+0.5)./(nOverall+1);

end