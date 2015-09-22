%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

%% plots the results produced by lassoFeatureSelectionWithReplication.m
%
% INPUTS_
%   resultfile_mat : matfile, containing the outputs of lassoFeatureSelectionWithReplication();
%   PRINTFLAG: boolean, indicating whether to export as pdfs
%   trueFeatures: a boolean vector indicating the true underlying features. these will be drawn as
%   solid, whereas all others are dashed!
function plotResults_LassoBootstrap_matfile(resultfile_mat,PRINTFLAG,trueFeatures)

load(resultfile_mat)

plotResults_LassoBootstrap(B,fitinfo,selectedFeatures,attributes ,PRINTFLAG,trueFeatures)




