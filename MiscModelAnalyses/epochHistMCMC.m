% Examine the firing rates in ephys vs whole cell recording data

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% ephys
nDatalist= 1;

figure;

barSeries   = 0:1:15;

for nData = nDatalist
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    load([TempDatDir DataSetList(nData).name '.mat'])
    m        = 2;
    
    dist        = 'Poisson';
    for nPlot   = 1:4
        nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
        subplot(m, m, nPlot)
        hold on;
        barData     = nan(length(nPeriodData), 1);
        for nUnit   = 1:length(nPeriodData)
            barData(nUnit) = mean(nPeriodData(nUnit).unit_yes_trial);
        end
        barSign     = 1;
        barHistWithDist(barData(:), dist, '', barSeries, 'b', barSign); 
        barData     = nan(length(nPeriodData), 1);
        for nUnit   = 1:length(nPeriodData)
            barData(nUnit) = mean(nPeriodData(nUnit).unit_no_trial);
        end
        barSign     = -1;
        barHistWithDist(barData(:), dist, '', barSeries, 'r', barSign); 
        hold off;
    end
end

load ([TempDatDir 'DataListC2SMCMCSingleTrialModel.mat']);
nDatalist= 1;
barSeries   = 0:1:15;
for nData = nDatalist
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    load([TempDatDir DataSetList(nData).name '.mat'])
    m        = 2;

    dist        = 'Poisson';
    for nPlot   = 1:4
        nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
        subplot(m, m, nPlot)
        hold on;
        barData     = nan(length(nPeriodData), 1);
        for nUnit   = 1:length(nPeriodData)
            barData(nUnit) = nanmean(nPeriodData(nUnit).unit_yes_trial)*params.frameRate;
        end
        barSign     = 1;
        barHistWithDist(barData(:), dist, '', barSeries, 'k', barSign); 
        barData     = nan(length(nPeriodData), 1);
        for nUnit   = 1:length(nPeriodData)
            barData(nUnit) = nanmean(nPeriodData(nUnit).unit_no_trial)*params.frameRate;
        end
        barSign     = -1;
        barHistWithDist(barData(:), dist, '', barSeries, 'g', barSign); 
        hold off;
    end
end

setPrint(8, 6, 'FR_distribution', 'pdf')