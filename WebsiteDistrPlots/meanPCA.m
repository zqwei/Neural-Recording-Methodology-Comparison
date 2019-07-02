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
ROCThres       = 0.50;

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
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
