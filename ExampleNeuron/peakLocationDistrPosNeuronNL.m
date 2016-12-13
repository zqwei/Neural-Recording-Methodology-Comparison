addpath('../Func');
setDir;

load ([TempDatDir 'DataListShuffle.mat']);
nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
spikeDataSet = nDataSet;   
params               = DataSetList(nData).params;

load('DeconvGPResults/S2CC2SMCMCSingleTrial_FineTunedModeled_GP43.mat', 'nDataSet')


numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
noProfileMatrix     = yesProfileMatrix;
positivePeak        = false(length(nDataSet));
for nUnit        = 1:length(nDataSet)
%     yesData      = mean(nDataSet(nUnit).decv_yes_trial, 1);
%     noData       = mean(nDataSet(nUnit).decv_no_trial, 1);
    yesData      = mean(nDataSet(nUnit).mcmc_yes_trial, 1);
    noData       = mean(nDataSet(nUnit).mcmc_no_trial, 1);
    maxData      = max([yesData, noData]);
    minData      = min([yesData, noData]);
    rData        = (maxData - minData);
    yesData      = (yesData - minData)/(maxData - minData);
    noData       = (noData - minData)/(maxData - minData);
    yesData(1:3) = nan;
    noData(1:3)  = nan;
    yesProfileMatrix(nUnit, :)    = yesData;
    noProfileMatrix(nUnit, :)     = noData;
    positivePeak(nUnit)       = mean(yesData(4:8)) <= mean(yesData(9:47)) ...
                               || mean(noData(4:8)) <= mean(noData(9:47));
end
actMat        = [yesProfileMatrix, noProfileMatrix];
actMat        = actMat(positivePeak, :);
[~, maxId]    = max(actMat, [], 2);

timeStep  = DataSetList(nData).params.timeSeries;
timeTag   = 8:60; % sample to response
numTime   = length(timeTag);
polein    = DataSetList(nData).params.polein;
poleout   = DataSetList(nData).params.poleout;

countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
% mean(countMaxId([8:77, 95:end]))
% disp(sqrt(mean((countMaxId([8:77, 95:end]) - 1/numTimeBin/2*100).^2)))
std(countMaxId([timeTag, timeTag+77]))%/mean(countMaxId([timeTag, timeTag+77]))
[bootstat,bootsam] = bootstrp(1000,@std,countMaxId([timeTag, timeTag+77]));
% mean(bootstat)
std(bootstat)
%     [h, p] = ttest(bootstat)
figure;
hold on;
%     stairs(timeStep, countMaxId(1:numTimeBin), '-', 'linewid', 1.0, 'color', [0.7 0 0])
bplot = bar(timeStep, countMaxId(1:numTimeBin), 1, 'facecolor', 'b', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
%     stairs(timeStep, countMaxId(1+numTimeBin:end), '-', 'linewid', 1.0, 'color', [0 0 0.7]);
bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 1, 'facecolor', 'r', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
xlim([timeStep(1) 2])
ylim([0 8])
gridxy ([polein, poleout, 0],[1/(numTimeBin-8)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
box off
ylabel('% Max peak')
xlabel('Time')
hold off
set(gca, 'TickDir', 'out')
setPrint(8, 3, 'SingleUnitsMaxLocationPosNeuron_S2C_C2S')

close all