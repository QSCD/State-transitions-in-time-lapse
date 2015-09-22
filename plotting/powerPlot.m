%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% plots the power as a function of sampleSize and Effect strength
% see Paper, Figure 5
%
% Inputs
%   powerMatrix: n x m x f Matrix of power values. 
%                n is the number of sampleSizes, 
%                m is the number of effect strength
%                f is the number of features
%   sampleSizes: 1 x n vector of different sample sizes used
%   effectStrength: 1 x m vector of different effect strength
%   featureNames: 1xf cell array of feature strings
%   subplotLayout:  (rownumbers, columnumbers)

function powerPlot(powerMatrix,sampleSizes,effectStrength,featureNames, subplotLayout )

nFeatures = size(powerMatrix,3);
assert(size(powerMatrix,1)==length(effectStrength))
assert(size(powerMatrix,2)==length(sampleSizes))

if ~exist('subplotLayout','var')
    rows  = floor(sqrt(nFeatures));
    cols = ceil(nFeatures/rows);
    subplotLayout = [rows,cols];
end


figure;

handles = subplot_no_margins(subplotLayout(1),subplotLayout(2));

[~, sortix_effect] = sort(effectStrength);
[~, sortix_sample] = sort(sampleSizes);

for i = 1:nFeatures
%     subplot(subplotLayout(1),subplotLayout(2),i)
    contourf(handles(i),sampleSizes(sortix_sample),log10(effectStrength(sortix_effect)),powerMatrix(sortix_effect,sortix_sample,i),0:0.2:1, 'Linestyle','None')
%     imagesc(powerArray(:,:,i))
    caxis(handles(i),[0 1])
    xlabel(handles(i),'Sample size')
    ylabel(handles(i),'log10(relative effect strength)')
    title(handles(i),featureNames{i});
    ylim(handles(i),[min(log10(effectStrength)) max(log10(effectStrength))])
    xlim(handles(i),[min(sampleSizes) max(sampleSizes)])
    
    set(gca,'XTick',[0 1000 2000 3000 4000])
    set(gca,'YTick',[-1,-0.5,0,0.5,1])
end

cm = gray;
colormap(cm(end:-1:1,:))
% colorbar

%get rid of the labels inside the plot
fixAxisLabels(handles)
end


function theAxes = subplot_no_margins(rows,cols)

baseAxis = gca;
pos = num2cell(get(baseAxis,'Position'));
[baseLeft,baseBottom, overallHeight, overallWidth] = deal(pos{:});

height = overallHeight/rows;
width = overallWidth/cols;

theAxes = [];

for i = 1:(rows*cols)
    [cRow,cCol] = ind2sub([rows cols], i);
    left = baseLeft + width*(cCol-1);
    bottom = baseBottom + height*(cRow-1);
   
   
    theAxes(i) = axes('Position',[left bottom width height]);  
    box on
end

delete(baseAxis)

%reshape such that we get i change first in rows

theAxes = reshape(theAxes,rows,cols);
%flip
theAxes =theAxes([end:-1:1],:);

%change from col first to row first
theAxes = theAxes';
theAxes = theAxes(:)';

%now the first entries is first row first col), second entry is first row
%second col

fixAxisLabels(theAxes)

end



function fixAxisLabels(theAxes)

%figure out what the lower and left border are
minLeft = Inf;
minBottom = Inf;
for i =1:length(theAxes)
   p = get(theAxes(i),'Position'); 
   
   l = p(1); b = p(2);
   
   minLeft = min(minLeft, l);
   minBottom = min(minBottom, b);
end

%for those plots at the border keep the axis, otherwise get rid of it
for i =1:length(theAxes)
   p = get(theAxes(i),'Position'); 
   
   l = p(1); b = p(2);   
   
   if l~=minLeft
        ylabel(theAxes(i),'')
        set(theAxes(i), 'YTickLabel', [],'YTick',[])       
   end

   if b~=minBottom
        xlabel(theAxes(i),'')
        set(theAxes(i), 'XTickLabel', [],'XTick',[])       
   end   
end
end

