%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% creates Fig5A,B: sample size vs inclusion probability
function sampleSizePlot(sampleSizes, means, trueFeaturesFlag)

% ColorSet = distinguishable_colors(14);
ColorSet = createCustomColorset();
StyleSet = {'-','x-','o-','d-'};

figure;hold on
set(gca, 'ColorOrder', ColorSet,'LineStyleOrder',StyleSet);

for i=1:size(means,2) %plot each feature iteratively
    if trueFeaturesFlag(i)==1
        style='-';
    else
        style='--';
    end
    
    % plot the density features with alternating markers
    markerStyle = 'None';
    markerSize = 5;
    if i>=5
       if mod(i,2)
            markerStyle = 'o';
       else
           markerStyle = 'None';
       end
    end
    plot(sampleSizes,means(:,i),'LineWidth',2, 'Color',ColorSet(i,:),'LineStyle',style, 'Marker', markerStyle, 'MarkerSize',markerSize, 'MarkerFaceColor',ColorSet(i,:))
end

set(gca,'XTick',sampleSizes)
set(gca,'YTick',[0,0.5,1])
xlabel('Sample size (onsets)')
ylabel('Mean inclusion probability')
hline(0.9,'r--')
xlim([min(sampleSizes),Inf])
legend({'t','x','y','cycle','$\phi_0$','$\phi_1$'  ,  '$\phi_{2}$'  ,  '$\phi_{3}$'  ,  '$\phi_{4}$' ,  '$\phi_{5}$'  ,  '$\phi_{6}$'  ,  '$\phi_{7}$','$\phi_{8}$','$\phi_{9}$'},'Interpreter','latex')
