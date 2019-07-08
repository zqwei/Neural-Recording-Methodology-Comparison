%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
xi = 0.0:0.01:1.0; 
color_ = [0.4667    0.6745    0.1882];

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 100;
numComps       = 3;
ROCThres       = 0.50;

fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);
    for nData      = 1:length(DataSetList)
        disp(DataSetList(nData).name);
        load([TempDatDir DataSetList(nData).name '.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
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
        for nPC = 1:3
            figure('visible', 'off');
            tmp_s = evMat(:, 2, nPC)./sum(evMat(:, :, nPC), 2);
            fs = ksdensity(tmp_s, xi);
            f_ = max(fs);
            plot(xi, fs/f_, '-', 'color', color_, 'linewid', 2);
            set(gca,'XTick',[0, 1], 'YTick', [0, 1], 'TickDir', 'out')
            ylabel('Prob. density')
            xlabel(['Frac. time content (PC' num2str(nPC) ')'])
            set(gca,'fontsize', 14)
            box off
            setPrint(8,3,[Result_ DataSetList(nData).name '_PC' num2str(nPC)], 'svg')
            close all;
        end
    end
end
