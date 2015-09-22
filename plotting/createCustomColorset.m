%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------
% creates the colormap used in the paper for the regularization plots
function ColorSet = createCustomColorset()
    nDensities = 10;
    %this is the red one
    black2red = [linspace(1,0,nDensities)' zeros(nDensities,1) zeros(nDensities,1)];
    otherColors = [0 0 1;0 1 0; 0 1 1; 1 0.8 0];
    ColorSet = [otherColors;black2red];

    %thats the green one
    black2green = [zeros(nDensities,1) linspace(1,0,nDensities)'  zeros(nDensities,1)];
    otherColors = [0 0 1;1 0 0; 0 1 1; 1 0.8 0];    
    ColorSet = [otherColors;black2green];   
    % save customColormap_red2black.mat ColorSet