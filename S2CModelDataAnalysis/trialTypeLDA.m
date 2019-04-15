%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
TempDatDir = '../Backups/TempDat_2019_01_28/';

numFold             = 30;
% load ([TempDatDir 'DataListS2CModel.mat']);
load([TempDatDir 'DataListS2CS1Model.mat'], 'DataSetList');
addNoise         = [1 0 0 0];

cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

numRandPickUnits    = 100;
numTrials           = numRandPickUnits*5;
numTestTrials       = numRandPickUnits*2;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

figure;

for nData             = 2%[1 3 4]
    if nData == 1
        load([TempDatDir 'Shuffle_Spikes.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    else
        load([TempDatDir DataSetList(nData).name '.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    end
    oldDataSet          = nDataSet;
    hold on
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    decodability          = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));
    for nFold    = 1:numFold
        trainingTargets     = [true(numTrainingTrials/2,1); false(numTrainingTrials/2,1)];
        trainingTargets     = trainingTargets(randperm(numTrainingTrials));
        testTargets         = [true(numTestTrials/2,1); false(numTestTrials/2,1)];
        testTargets         = testTargets(randperm(numTestTrials));
        totTargets          = [testTargets; trainingTargets];

        trainingDecisions   = trainingTargets(randperm(numTrainingTrials));
        testDecisions       = testTargets(randperm(numTestTrials));
        totDecisions        = [testDecisions; trainingDecisions];

        randPickUnits       = randperm(length(nDataSet));
        randPickUnits       = randPickUnits(1:numRandPickUnits);

        nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
        decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
    end

    meandecodability = mean(decodability,1);
%     meanValue        = mean(meandecodability(1:8));
%     meandecodability = (meandecodability - meanValue)/(1-meanValue)*0.5+0.5;

%     shadedErrorBar(DataSetList(nData).params.timeSeries, meandecodability,...
%         std(decodability, 1)/sqrt(numFold),...
%         {'-', 'linewid', 1.0, 'color', cmap(nData,:)}, 0.5);
    plot(DataSetList(nData).params.timeSeries, meandecodability, 'k');
    plot(DataSetList(nData).params.timeSeries, meandecodability+std(decodability, 1)/sqrt(numFold), 'k');
    plot(DataSetList(nData).params.timeSeries, meandecodability-std(decodability, 1)/sqrt(numFold), 'k');  
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    gridxy ([DataSetList(nData).params.polein, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    set(gca, 'TickDir', 'out')
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
end

setPrint(8, 6, [PlotDir 'S2CModel/CollectedUnitsDecodabilityROC_Summary'], 'pdf')
margNames = {'Spike', '', 'S2C long decay', 'S2C short decay'};

figure;
hold on
for nColor = [1 3 4]
    plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
    text(1, nColor, margNames{nColor})
end
xlim([0 10])
hold off
axis off
setPrint(3, 2, [PlotDir 'S2CModel/CollectedUnitsDecodabilityROC_SummaryLabel'])

close all;
