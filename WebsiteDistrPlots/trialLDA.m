%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
file_ext = {'Decode_Sample', 'Decode_Delay'};
x_labels = {'Sample-epoch decodability', 'Delay-epoch decodability'};
xi = 0.5:0.01:1.0; 
color_ = [0.7176    0.2745    1.0000];

numTrials      = 100;
numRandPickUnits = 50;
numTrainingTrials = 300;
numTestTrials     = 300;
numFold        = 100;
ROCThres       = 0.50;

fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);

    for nData      = 1:length(DataSetList)
        disp(DataSetList(nData).name);
        load([TempDatDir DataSetList(nData).name '.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
        decodability = zeros(numFold, size(nDataSet(nData).unit_yes_trial,2));
        params      = DataSetList(nData).params;
        timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
        valid_cell = true(length(nDataSet), 1);
        for nNeuron = 1:length(nDataSet)
            if any(isnan(nDataSet(nNeuron).unit_yes_trial), 'all') || any(isnan(nDataSet(nNeuron).unit_no_trial), 'all')
                valid_cell(nNeuron) = false;
            end
        end
        nDataSet = nDataSet(valid_cell);
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
            if numRandPickUnits < length(nDataSet)
                randPickUnits   = randPickUnits(1:numRandPickUnits);
            end
            if length(nDataSet) < 10
                continue;
            end
            nSessionData        = shuffleSessionData(nDataSet(randPickUnits), totTargets, numTestTrials);
            decodability(nFold,:) = decodabilityLDA(nSessionData, trainingTargets, testTargets);
        end
        for nPC = 1:2
            figure('visible', 'off');
            tmp_s = mean(decodability(:, timePoints(nPC+1):timePoints(nPC+2)), 2);
            fs = ksdensity(tmp_s, xi);
            f_ = max(fs);
            plot(xi, fs/f_, '-', 'color', color_, 'linewid', 2);
            xlim([0.5, 1])
            ylim([0, 1])
            set(gca,'XTick',[0.5, 1], 'YTick', [0, 1], 'TickDir', 'out')
            ylabel('Prob. density')
            xlabel(x_labels{nPC})
            set(gca,'fontsize', 14)
            box off
            setPrint(8,4.5,[Result_ DataSetList(nData).name '_' file_ext{nPC}], 'svg')
            close all;
        end
    end
end
