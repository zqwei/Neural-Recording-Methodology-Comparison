addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

timeTag       = 8;

for nData     = [3 4]
    load(['S2CC2S_' DataSetList(nData).name '.mat'])  
    spkDataSet = nDataSet;
    for nCell  = 1:length(nDataSet)
        spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).mcmc_yes_trial;
        spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).mcmc_no_trial;
    end
    nDataSet            = spkDataSet;
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

    timeStep  = DataSetList(nData).params.timeSeries;
    numTime   = length(timeStep);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    
    
    countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
%     std(countMaxId([timeTag, timeTag+77]))%/mean(countMaxId([timeTag, timeTag+77]))
    [bootstat,bootsam] = bootstrp(1000,@mean,countMaxId([timeTag, timeTag+77]));
    mean(bootstat)
    std(bootstat)
%     [h, p] = ttest(bootstat)
    figure;
    hold on;
    bplot = bar(timeStep, countMaxId(1:numTimeBin), 'facecolor', 'b', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
    bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 'facecolor', 'r', 'edgecolor', 'none');
    bplot.FaceAlpha = 0.5;
    xlim([timeStep(1) 2])
    ylim([0 8])
    gridxy ([polein, poleout, 0],[1/(numTimeBin-8)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
    box off
    ylabel('% Max peak')
    xlabel('Time')
    hold off
    setPrint(8, 3, [PlotDir 'S2CModel\SingleUnitsMaxLocationPosNeuronS2CC2S_' DataSetList(nData).name])
end

% close all


fName = 'FineTunedModeled_GP517';
load(['S2CC2S_' fName '.mat'])  
spkDataSet = nDataSet;
for nCell  = 1:length(nDataSet)
    spkDataSet(nCell).unit_yes_trial = spkDataSet(nCell).mcmc_yes_trial;
    spkDataSet(nCell).unit_no_trial  = spkDataSet(nCell).mcmc_no_trial;
end
nDataSet            = spkDataSet;
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

timeStep  = DataSetList(nData).params.timeSeries;
numTime   = length(timeStep);
polein    = DataSetList(nData).params.polein;
poleout   = DataSetList(nData).params.poleout;


countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
%     std(countMaxId([timeTag, timeTag+77]))%/mean(countMaxId([timeTag, timeTag+77]))
[bootstat,bootsam] = bootstrp(1000,@mean,countMaxId([timeTag, timeTag+77]));
mean(bootstat)
std(bootstat)
%     [h, p] = ttest(bootstat)
figure;
hold on;
bplot = bar(timeStep, countMaxId(1:numTimeBin), 'facecolor', 'b', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 'facecolor', 'r', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
xlim([timeStep(1) 2])
ylim([0 8])
gridxy ([polein, poleout, 0],[1/(numTimeBin-8)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
box off
ylabel('% Max peak')
xlabel('Time')
hold off
setPrint(8, 3, [PlotDir 'S2CModel\SingleUnitsMaxLocationPosNeuronS2CC2S_' fName])