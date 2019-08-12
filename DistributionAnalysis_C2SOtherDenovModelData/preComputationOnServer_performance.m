addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');

numComps       = 3;
combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
numTrials      = 100;
numFold        = 1;
numTrialsLDA   = 500;
numTestTrials  = 200;
numTrainingTrials = numTrialsLDA - numTestTrials;
numRandPickUnits  = 50;
ROCThres          = ones(10,1)*0.5;
ROCThres(10)      = 0.7;

compiled_results  = [];
for nData         = [1 10]  
    load([TempDatDir DataSetList(nData).name '.mat']);   
    analysisMat   = repmat(struct('nParaSet',1, 'pCa0', 1, ...
                                'pTaud', 1, 'sizeGroup', 1),100, 1);
    yesData       = nDataSet(nData).unit_no_trial;
    numUnit       = length(nDataSet);
    numT          = size(yesData, 2);
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    selectedNeuronalIndex = selectedHighROCneurons(nDataSet, DataSetList(nData).params, ROCThres(nData), selectedNeuronalIndex);
    depth_list          = [nDataSet.depth_in_um]';
    selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;
    spikeDataSet  = nDataSet(selectedNeuronalIndex);
    numUnit       = length(spikeDataSet);

    for nParaSet  = 1:1000              
        disp(nParaSet)        
        nDataSet  = spikeDataSet(rand(numUnit, 1)<0.9);        
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
        analysisMat(nParaSet).PCAVar       = PCAVar;
        analysisMat(nParaSet).decodability = decodability;
    end
    tmp.name         = DataSetList(nData).name;
    tmp.analysisMat  = analysisMat;
    compiled_results = [compiled_results, tmp];
end


save('Results_compiled_raw_6f.mat', 'compiled_results');