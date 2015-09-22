%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

% figure out which selected features correspond to density features (\phi_k)!
function kernelWeights = extractDensityFeatures(selectedFeatureWeights,attributes,attributes_selected,NOCONSTANT_FLAG)

densFeatures_names = attributes(5:end);

if NOCONSTANT_FLAG == 0
    featureWeights = selectedFeatureWeights(2:end);
else
    featureWeights = selectedFeatureWeights(1:end);
end
kernelWeights = helper(featureWeights,densFeatures_names,attributes_selected);
end


function kernelWeights = helper(featureWeights,densFeatures_names,selectedFeaturesCellArray)

    kernelWeights=  zeros(size(densFeatures_names));
    %which of those are density related
    for i = 1:length(densFeatures_names)
        ix_tmp = strcmp(selectedFeaturesCellArray,densFeatures_names{i}); %is there a selected feature that corresponds to denstiy?
        if any(ix_tmp)
            kernelWeights(i)= featureWeights(ix_tmp); % if yes, select its coefficient
        else
            kernelWeights(i)=0;
        end
        
    end
end