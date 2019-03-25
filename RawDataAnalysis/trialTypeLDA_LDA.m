%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
ROCThres            = 0.50;
cmap                = cbrewer('div', 'Spectral', 128, 'cubic');
numFold             = 30;

nDataList      = [1 3 4 13 10 11 12];

for mData      = 1:length(nDataList)
    nData      = nDataList(mData);
    if nData == 13
        load([TempDatDir DataSetList(1).name '.mat'])
        depth_list = [nDataSet.depth_in_um];
        nDataSet = nDataSet(depth_list < 471);
        nData = 1;
        nName  = 'Ephys_6f';
    else
        load([TempDatDir DataSetList(nData).name '.mat'])
        nName  = DataSetList(nData).name;
    end
    
    params                = DataSetList(nData).params;
    numT                  = length(params.timeSeries);
    numRandPickUnits      = length(nDataSet);
    numTrials             = numRandPickUnits*3;
    totTargets            = [true(numTrials,1); false(numTrials,1)];
    corrMat               = zeros(numFold, numT, numT);

    for nFold             = 1:numFold
        nSessionData          = shuffleSessionData(nDataSet, totTargets, numTrials*2);
        nSessionData          = normalizationDim(nSessionData, 2);
        coeffs                = coeffLDA(nSessionData, totTargets);
        corrMat(nFold, :, :)  = coeffs'*coeffs;
    end

    figure;
    imagesc(params.timeSeries, params.timeSeries, squeeze(mean(corrMat, 1)));
    xlim([min(params.timeSeries) max(params.timeSeries)]);
    ylim([min(params.timeSeries) max(params.timeSeries)]);
    caxis([0 1]);
    axis xy;
    gridxy ([params.polein, params.poleout, 0],[params.polein, params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    colormap(cmap)
    xlabel('LDA Time (s)')
    ylabel('LDA Time (s)')
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, ['SimilarityLDALDA_' nName], 'pdf')
end


load ([TempDatDir 'DataListS2CShuffle.mat']);
nDataList      = [3 6 15 18];

for mData      = 1:length(nDataList)
    nData      = nDataList(mData);
    load([TempDatDir DataSetList(nData).name '.mat'])
    nName  = DataSetList(nData).name;
    
    params                = DataSetList(nData).params;
    numT                  = length(params.timeSeries);
    numRandPickUnits      = length(nDataSet);
    numTrials             = numRandPickUnits*3;
    totTargets            = [true(numTrials,1); false(numTrials,1)];
    corrMat               = zeros(numFold, numT, numT);

    for nFold             = 1:numFold
        nSessionData          = shuffleSessionData(nDataSet, totTargets, numTrials*2);
        nSessionData          = normalizationDim(nSessionData, 2);
        coeffs                = coeffLDA(nSessionData, totTargets);
        corrMat(nFold, :, :)  = coeffs'*coeffs;
    end

    figure;
    imagesc(params.timeSeries, params.timeSeries, squeeze(mean(corrMat, 1)));
    xlim([min(params.timeSeries) max(params.timeSeries)]);
    ylim([min(params.timeSeries) max(params.timeSeries)]);
    caxis([0 1]);
    axis xy;
    gridxy ([params.polein, params.poleout, 0],[params.polein, params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
    box off;
    hold off;
    colormap(cmap)
    xlabel('LDA Time (s)')
    ylabel('LDA Time (s)')
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, ['SimilarityLDALDA_' nName], 'pdf')
end