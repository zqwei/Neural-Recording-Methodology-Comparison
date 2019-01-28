%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Same ROC
% Different number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('../Func');
setDir;

numFold             = 10;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
ROCThres            = 0.5;

load ([TempDatDir 'DataListShuffle.mat']);
addNoise            = zeros(length(DataSetList), 1);

cmap                = cbrewer('div', 'Spectral', 128, 'cubic');


if ~exist([PlotDir '/CollectedUnitsDecodability'],'dir')
    mkdir([PlotDir '/CollectedUnitsDecodability'])
end

stepSize            = 50; % 15 for the last data.

for nData      = 1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    oldDataSet               = nDataSet;

    for numRandPickUnits      = 1:length(stepSize)        
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex(~neuronRemoveList)';% & [DataSetList(nData).cellinfo(:).cellType] == 1;
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
            num_                = min(stepSize(numRandPickUnits), length(randPickUnits)-1);
            randPickUnits       = randPickUnits(1:num_);

            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData +randn(size(nSessionData))*1e-3/sqrt(numTrials)* addNoise(nData), trainingTargets, testTargets);
        end
        figure;
        hold on
        shadedErrorBar(DataSetList(nData).params.timeSeries, mean(decodability,1),...
            std(decodability, 1)/sqrt(numFold),...
            {'-', 'linewid', 1.0, 'color', 'k'}, 0.5);
        plot(DataSetList(nData).params.timeSeries, mean(decodability), '-k', 'linewid', 2.0);
        xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
        ylim([0.5 1])
        gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
        box off;
        hold off;
        xlabel('Time (s)');
        ylabel('# units');
        set(gca, 'TickDir', 'out')
        setPrint(8, 6, [PlotDir 'WebsitePlots/' DataSetList(nData).name '_decodability'], 'svg')
    end    
end

close all