%--------------------------------------------------
%   Author: Michael Strasser
%   Helmholtz Zentrum München
%   Institute of Computational Biology
%   http://www.helmholtz-muenchen.de/icb/
%   09/2015
%--------------------------------------------------

function createFeatureMatrix(filename_trees_notPruned_withDensities, output_name, markerField)

disp('loading')
load (filename_trees_notPruned_withDensities)

disp('pruning')
%preprocess: prune trees at differentiation event, add the \phi_k, add cell cycle
pC_time = @(t)t.(markerField)==1;
g = @(mT)tUtil_pruneMichiTree_TIME(mT,pC_time);
trees = cellfun(g, trees_notPruned_withDensities,'UniformOutput',0);

disp('adding rings/cctime')
trees = add_RingDensityFeatures(trees);
trees = addTimeFromCCstart(trees);

attributes = {'timepoint','cycleTime','X','Y','ring_density_0','ring_density_10','ring_density_20','ring_density_30','ring_density_40','ring_density_50','ring_density_60','ring_density_70','ring_density_80','ring_density_90'};


%% gather the data/labels from the trees
%to prealloc the arrays, figure out how many onsets, how many non onsets
% htis has to be done on the pruned trees!@!
onsets = (cellfun(@(t)sum(t.(markerField)==1), trees));
noOnsets = cellfun(@(t)sum(t.(markerField)==0), trees);
X_prealloc = zeros(sum(onsets), length(attributes) );
X_not_prealloc = zeros(sum(noOnsets), length(attributes) );
counter_did = 1; %next empty element in the preloc array
counter_didnot = 1;


for i = 1:length(trees)
    fprintf('creating feature matrix %d/%d\n', i, length(trees))
    cT = trees{i};
    ix_did = cT.(markerField)==1;

    for j=1:length(attributes)
        cA= attributes{j};
        dataDid = cT.(cA)(ix_did)';
        dataDidnot = cT.(cA)(~ix_did)';

        
        ixEntries_did = counter_did:(counter_did+length(dataDid)-1);
        X_prealloc(ixEntries_did,j) = dataDid;
        
        ixEntries_didnot = counter_didnot:(counter_didnot+length(dataDidnot)-1);
        X_not_prealloc(ixEntries_didnot,j) = dataDidnot;
    end
        counter_did = counter_did+length(dataDid);
        counter_didnot = counter_didnot+length(dataDidnot);       
end

%rename the preallocs into the old names
X = X_prealloc;
X_not = X_not_prealloc;
treeId_cellId_array = nan; %TODO reimplement this with prealloc as well
treeId_cellId_array_not = nan; %TODO reimplement this with prealloc as well

%remove rows with nans
ix1 = any(isnan(X),2);
X = X(~ix1,:);


ix2 = any(isnan(X_not),2);
X_not = X_not(~ix2,:);


disp('saving')
save([output_name '.mat']', 'X', 'X_not', 'attributes', 'treeId_cellId_array_not', 'treeId_cellId_array') 
