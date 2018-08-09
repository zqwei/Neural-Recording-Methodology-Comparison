%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numComps       = 3;

if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end

ROCThres            = 0.65;


for nData              = 10%[1 3 4 10]
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

    
    params      = DataSetList(nData).params;
    numTime     = length(params.timeSeries);
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaFiringRatesAverage = zeros(numComps, 2, numTime);
    firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
    [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
    pcaFiringRatesAverage(:, 1, :) = score(1:numTime, :)';
    pcaFiringRatesAverage(:, 2, :) = score(numTime+1:end, :)';
    
    figure;
    for nPlot           = 1:numComps
        subplot(1, numComps, nPlot)
        hold on
        plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 1, :)), '-b', 'linewid', 2);
        plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 2, :)), '-r', 'linewid', 2);
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        hold off
        box off
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        xlabel('Time (s)')
        ylabel(['PC' num2str(nPlot) ' score'])  
        set(gca, 'TickDir', 'out')
    end
    
    setPrint(8*3, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCATrace_' DataSetList(nData).name])
    
end

% ROCThres = 0.55;
% % different ROC
% for nData              = [1 3 4]
%     load([TempDatDir DataSetList(nData).name '.mat']);
%     selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
%     selectedNeuronalIndex = selectedHighLocalROCneurons(nDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
%     nDataSet              = nDataSet(selectedNeuronalIndex);
%     firingRates        = generateDPCAData(nDataSet, numTrials);
%     firingRatesAverage = nanmean(firingRates, ndims(firingRates));
%     pcaFiringRatesAverage = zeros(numComps, 2, 77);
%     firingRatesAverage = [squeeze(firingRatesAverage(:, 1, :)), squeeze(firingRatesAverage(:, 2, :));];
%     [~,score,~]        = pca(firingRatesAverage', 'NumComponents', numComps);
%     pcaFiringRatesAverage(:, 1, :) = score(1:77, :)';
%     pcaFiringRatesAverage(:, 2, :) = score(78:end, :)';
%     
%     figure;
%     for nPlot           = 1:numComps
%         subplot(1, numComps, nPlot)
%         hold on
%         plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 1, :)), '-r', 'linewid', 2);
%         plot(DataSetList(nData).params.timeSeries, squeeze(pcaFiringRatesAverage(nPlot, 2, :)), '-b', 'linewid', 2);
%         gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
%         hold off
%         box off
%         xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
%         xlabel('Time (s)')
%         ylabel(['PC' num2str(nPlot) ' score'])  
%     end
%     
%     setPrint(8*3, 6, [PlotDir 'CollectedUnitsPCA/CollectedUnitsPCALOCThresTrace_' DataSetList(nData).name])
%     
% end
% 
% close all