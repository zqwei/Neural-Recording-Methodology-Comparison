%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.  Collected population decision decodability over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V1_6_1
% The same amount of units in analysis




addpath('../Func');
setDir;

numFold             = 50;
numTrials           = 500;
numTestTrials       = 200;
numTrainingTrials   = numTrials - numTestTrials;
numRandPickUnits    = 100;

load ([TempDatDir 'DataListModeled.mat']);
addNoise         = [1 0 0 0 0 0];

if ~exist([PlotDir '/Collected_ModelUnits_Decodability'],'dir')
    mkdir([PlotDir '/Collected_ModelUnits_Decodability'])
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    figure;
    nDataSet     = nDataSet(selectedNeuronalIndex);
    decodability = zeros(numFold, size(nDataSet(1).unit_yes_trial,2));
    
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
    hold on
    plot(DataSetList(nData).params.timeSeries, mean(decodability,1), 'k', 'linewid',1);
    xlim([min(DataSetList(nData).params.timeSeries) max(DataSetList(nData).params.timeSeries)]);
    ylim([0.5 1])
    gridxy ([DataSetList(nData).params.polein, DataSetList(nData).params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 0.5)
    box off;
    hold off;
    xlabel('Time (s)');
    ylabel('Decodability');
    setPrint(4, 3, [PlotDir 'Collected_ModelUnits_Decodability/Collected_Units_Decodability_FixedNumberUnits_' DataSetList(nData).name], 'pdf')
end
close all;