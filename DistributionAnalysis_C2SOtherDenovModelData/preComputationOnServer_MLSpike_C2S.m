addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
dataSetNames   = {'ModelSpikeMLSpike_Deconv_Ca_Fast_SShort_Delay','ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay_Virus', 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay'};
indexDatasets  = [10, 4, 3];

numComps       = 3;
combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 1;
numTrialsLDA   = 500;
numTestTrials  = 200;
numTrainingTrials = numTrialsLDA - numTestTrials;
numRandPickUnits  = 50;

for nData         = 1:length(dataSetNames)        
    load([TempDatDir dataSetNames{nData} '.mat']);   
    analysisMat   = repmat(struct('nParaSet',1, 'pCa0', 1, ...
                                'pTaud', 1, 'sizeGroup', 1),100, 1);
    yesData       = nDataSet(1).unit_no_trial;
    numUnit       = length(nDataSet);
    numT          = size(yesData, 2);
    timeTag       = 8:(numT-3);
    spikeDataSet  = nDataSet;
    params        = DataSetList(indexDatasets(nData)).params;
    
    for nParaSet  = 1:100              
        disp(nParaSet)        
        nDataSet  = spikeDataSet(rand(numUnit, 1)<0.9);        
        % cell type data
        unitGroup = getLogPValueTscoreSpikeTime(nDataSet, params);       
        sizeGroup = histcounts(unitGroup, 0:3); % record this for cell type data
        
        % peakiness
        numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
        yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
        noProfileMatrix     = yesProfileMatrix;
        positivePeak        = false(length(nDataSet));
        for nUnit        = 1:length(nDataSet)
            yesData      = mean(nDataSet(nUnit).unit_yes_trial);
            noData       = mean(nDataSet(nUnit).unit_no_trial);
            maxData      = max([yesData, noData]);
            minData      = min([yesData, noData]);
            rData        = (maxData - minData);
            yesData      = (yesData - minData)/(maxData - minData);
            noData       = (noData - minData)/(maxData - minData);
            yesProfileMatrix(nUnit, :)    = yesData;
            noProfileMatrix(nUnit, :)     = noData;
            positivePeak(nUnit)       = mean(yesData(1:8)) <= mean(yesData(9:47)) ...
                                       || mean(noData(1:8)) <= mean(noData(9:47));
        end
        actMat        = [yesProfileMatrix, noProfileMatrix];
        actMat        = actMat(positivePeak, :);
        [~, maxId]    = max(actMat, [], 2);
        countMaxId    = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
        peakiness     = std(countMaxId([timeTag, timeTag+numT]));
        
        % pca
        evMat              = zeros(numFold, length(combinedParams), numComps);
        firingRates        = generateDPCAData(nDataSet, numTrials);
        firingRatesAverage = nanmean(firingRates, ndims(firingRates));
        pcaX               = firingRatesAverage(:,:);
        firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
        pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
        Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
        totalVar           = sum(sum(pcaX.^2));
        [~, ~, Wpca] = svd(pcaX');
        PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
        for i=1:length(Xmargs)
            PCAmargVar(i,:)= sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
        end
        PCAVar             = PCAmargVar(:, 1:numComps)';
        
        %lda
        decodability            = nan(numFold, size(nDataSet(1).unit_yes_trial,2));        
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
        analysisMat(nParaSet).nParaSet     = nParaSet;
        analysisMat(nParaSet).sizeGroup    = sizeGroup;
        analysisMat(nParaSet).peakiness    = peakiness;
        analysisMat(nParaSet).PCAVar       = PCAVar;
        analysisMat(nParaSet).decodability = decodability;
    end
    
    save([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNames{nData} '.mat'], 'analysisMat');
end