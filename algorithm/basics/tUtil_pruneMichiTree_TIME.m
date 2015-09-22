%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% prunes the tree exactly at the timepoints where the condition is fullfilled the first time
% (the timepoint where the condition is fullfilled the first time stays in the tree)
function [prunedTree, leaves_ofPruned ]= tUtil_pruneMichiTree_TIME(currentTree,pruningCond, startingCell_dirtyHack)

   leaves_ofPruned = [];
   
   if ~exist('startingCell_dirtyHack','var')
        nodes = [1];
        assert(any(currentTree.cellNr==1),'pruning function doesnt work for cells without a root')
   else
       st=dbstack;
       if ~strcmp(st(2).name, 'spliceTree') %this functionality is only supposed to be used from within this spliceTree function, if called from another function issue this warning
        warning('this is some dirty hack to make pruning work for trees that dont have root=1. not really tested, only used in treeSplicing')
       end
       nodes = startingCell_dirtyHack;
   end
   
   retainedNodes = [];
   
   
   retained_indices = zeros(size(currentTree.cellNr)); %which datapoints are retained after pruning

   pC = pruningCond(currentTree);
   assert(length(pC)==length(currentTree.cellNr)); % the pruning condition must return one value per timepoint for a given cell
   
   allCells = unique(currentTree.cellNr);

   while (~isempty(nodes))
       nextLevel = [];
       for i=1:length(nodes)
           currentNode = nodes(i);
           
           assert( currentNode ~= 0)
           
           retainedNodes= [retainedNodes currentNode];
           
           cellix = currentTree.cellNr == currentNode;           
           
           pruneThisCell = any(pC & cellix); %pruning cond fulfilled in this cell?
           % if not, process children
           if ~pruneThisCell 
               
               % all timepoints of the current cell are retained
               retained_indices = retained_indices+cellix;
               
               %left
               if any(allCells == 2.*currentNode)
                   nextLevel = [nextLevel 2.*currentNode ];
               end
               %right
               if any(allCells== 2.*currentNode+1)
                   nextLevel = [nextLevel 2.*currentNode+1 ];
               end
               
           else
               % only those datapoints are retained up to (and including) the time when the
               % condition is fulfilled
               firstPoint_ofCell = find(cellix,1,'first');
               first_pC_point = find(pC & cellix,1,'first');
               retained_indices(firstPoint_ofCell:first_pC_point)=1;
               
               %remember the leave
               leaves_ofPruned = [leaves_ofPruned currentNode];
           end
       end
          
       nodes = nextLevel;
   end
   assert(all(retained_indices<2))
   
   %% put together the pruned tree from the retained indices   
   prunedTree = tUtil_treeIndexing(currentTree,retained_indices==1);
end