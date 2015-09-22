%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% does the barplot of fig6A to compare the observed and expected sister
% freqs
function plot_sister_correlation_figure(O_true,E_true,E_impaired, pval_true,pval_impaired  )

figure;
colormap([0 0 0;0.7 0.7 0.7; 1 1 1])
bar([O_true;E_true; E_impaired]')
set(gca,'XTickLabel',{'both','one','none'})
legend({'Observed','Expected','Expected impaired'})
box off
ylabel('Absolute frequency of sister subtrees')
xlabel('Sister subtrees')
title(['pval expect ' num2str(pval_true) ' pval impaired ' num2str(pval_impaired)])
set(gca,'yscale','log')
ylim([100, 2e4])
box off

end