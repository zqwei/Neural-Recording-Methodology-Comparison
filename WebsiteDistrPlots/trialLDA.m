%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
% load ([TempDatDir 'DataListC2SShuffle.mat']);
% load ([TempDatDir 'DataListS2CShuffle.mat']);
file_ext = {'Peakiness'};
xi = 0.0:0.1:2.0; 
color_ = [0.4313    0.4313    0.4313];

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 100;

numComps       = 10;
cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

ROCThres            = 0.50;


for nData      = [1 3 4]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';
    end
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2)); 
    
    evMat              = zeros(numFold, length(combinedParams), numComps);
    for nFold          = 1:numFold
        firingRates        = generateDPCAData(nDataSet, numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));
        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        totalVar           = sum(sum(pcaX.^2));
        [~, S, Wpca] = svd(pcaX');

        PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
        for i=1:length(Xmargs)
            PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
        end
        evMat(nFold, :, :) = PCAmargVar(:, 1:numComps);
    end
    if nData == 1
        refEphys.pca = evMat;
        disp(mean(evMat(:, :, 1)))
    end
    if nData == 3
        performanceMat(2).pca = evMat;
        disp(mean(evMat(:, :, 1)))
    end
    if nData == 4
        performanceMat(1).pca = evMat;
        disp(mean(evMat(:, :, 1)))
    end
end
