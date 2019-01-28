addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

for nData     = 1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    
    disp(DataSetList(nData).name)
    params      = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    
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
        positivePeak(nUnit)       = mean(yesData(timePoints(1):timePoints(2))) <= mean(yesData(timePoints(2):timePoints(4))) ...
                                   || mean(noData(timePoints(1):timePoints(2))) <= mean(noData(timePoints(2):timePoints(4)));

    end
    actMat        = [yesProfileMatrix, noProfileMatrix];
    actMat        = actMat(positivePeak, :);
    [~, maxId]    = max(actMat, [], 2);

    timeStep  = DataSetList(nData).params.timeSeries;
    timeTag   = timePoints(2):timePoints(4)+13;
    numTime   = length(timeTag);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
    figure;
    hold on;
    bplot = bar(timeStep, countMaxId(1:numTimeBin), 1, 'facecolor', 'b', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
    bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 1, 'facecolor', 'r', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
    xlim([timeStep(1) params.timeSeries(end)]);
    ylim([0 13])
    preBin = sum(DataSetList(nData).params.timeSeries<polein);
    gridxy ([polein, poleout, 0],[1/(numTimeBin-preBin)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
    box off
    set(gca,'TickDir','out')
    ylabel('% Max peak')
    xlabel('Time')
    hold off
    setPrint(8, 3, [PlotDir 'WebsitePlots\' DataSetList(nData).name '_peakness'], 'svg')
end

close all