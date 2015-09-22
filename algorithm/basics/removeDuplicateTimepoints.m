%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% depending on the simulation engine, daughters might start at the same timepoint as the mother
% ended -> this timepoint is then counted to often (for stat. analysis); remove them
function newTree = removeDuplicateTimepoints(tree)

    removeIx= zeros(1,length(tree.cellNr));
    
    uniqueCells = unique(tree.cellNr);
    
    for i=1:length(uniqueCells)
        cC= uniqueCells(i);
        
        %check if the last timepoint of this cell is present in daughters
        %of this cell
        lastTP = max(tree.absoluteTime(tree.cellNr==cC));
        d1 = cC.*2; d2 = cC.*2+1;
        if any(uniqueCells==d1)
            firstTP = min(tree.absoluteTime(tree.cellNr==d1));
            if firstTP==lastTP;
                removeIx(tree.cellNr==d1 & tree.absoluteTime==firstTP) =1; %mark to remove this one
            end
        end
        
        if any(uniqueCells==d2)
            firstTP = min(tree.absoluteTime(tree.cellNr==d2));
            if firstTP==lastTP;
                removeIx(tree.cellNr==d2 & tree.absoluteTime==firstTP) =1; %mark to remove this one
            end
        end        
    end
    
    %take those not marked for removal
    newTree = tUtil_treeIndexing(tree,removeIx==0);

end