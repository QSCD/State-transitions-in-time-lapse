%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% given weights of the \phi_k
% plot the resconstructed kernel
function reconstruct_Densitykernel_from_rings(weights) 

    reconstr_kernel = zeros(1,105);
    ringLocation = [0 10 20 30 40 50 60 70 80 90]; %the position of the inner part of the ring
    dR = 10;
    
    assert(length(weights) == length(ringLocation))
    
    for i = 1:length(weights)
        
        currentRingLoc = ringLocation(i) +1 ;
        %the function is only non  0 if inside [ringlocation,ringlocation+dR]
        currentFunction = zeros(1,105);
        currentFunction(currentRingLoc:(currentRingLoc+dR-1))=1; 
        reconstr_kernel = reconstr_kernel + weights(i).*currentFunction;
    end
       
    figure;hold on
    bar(ringLocation,abs(weights),'Style','histc') %,'LineStyle','none'
    bar(-(ringLocation+10),abs(weights),'Style','histc') %histc left aligns so we have to shift by 10
    xlabel('Distance')
    ylabel('Density weight')
end