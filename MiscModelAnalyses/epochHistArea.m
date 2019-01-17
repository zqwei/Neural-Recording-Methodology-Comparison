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
    timePoints  = [timePoints(2), timePoints(end)];
    load([TempDatDir DataSetList(nData).name '.mat'])
    m        = 2;
    
    dist        = 'Poisson';
    nPlot       = 1;
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
    hold on;
    barData     = nan(length(nPeriodData), 2);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit,1) = mean(nPeriodData(nUnit).unit_yes_trial);
    end
    barSign     = 1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'b', barSign); 
%     barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 2) = mean(nPeriodData(nUnit).unit_no_trial);
    end
%     barSign     = -1;
    barHistWithDist(mean(barData, 2), dist, '', barSeries, 'r', barSign); 
    hold off;
end

% ca_factor = [9.5005, 1.9279, 5.4207];
% 
% ca_factor = ca_factor./[60, 55, 55] * 15;

ca_factor = [2.3701, 0.3649, 1.2353];

hold on;
for nData = 4
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    timePoints  = [timePoints(2), timePoints(end)];
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    m        = 2;
    
    dist        = 'Poisson';
    nPlot       = 1;
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
    hold on;
    barData     = nan(length(nPeriodData), 2);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 1) = sum(max(nPeriodData(nUnit).unit_yes_trial, 0))/ca_factor(1)/(params.timeSeries(end)-params.polein);
    end
%     barSign     = 1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'g', barSign); 
%     barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 2) = sum(max(nPeriodData(nUnit).unit_no_trial, 0))/ca_factor(1)/(params.timeSeries(end)-params.polein);
    end
%     barSign     = -1;
    barHistWithDist(mean(barData, 2), dist, '', barSeries, 'k', barSign); 
    hold off;
end

% nDatalist= 1;
% barSeries   = 0:1:15;
% for nData = nDatalist
%     params   = DataSetList(nData).params;
%     timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
%     load([TempDatDir 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay.mat'])
%     m        = 2;
% 
%     dist        = 'Poisson';
%     for nPlot   = 1:4
%         nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
%         subplot(m, m, nPlot)
%         hold on;
%         barData     = nan(length(nPeriodData), 1);
%         for nUnit   = 1:length(nPeriodData)
%             barData(nUnit) = nanmean(nPeriodData(nUnit).unit_yes_trial)*params.frameRate;
%         end
%         barSign     = 1;
%         barHistWithDist(barData(:), dist, '', barSeries, 'k', barSign); 
%         barData     = nan(length(nPeriodData), 1);
%         for nUnit   = 1:length(nPeriodData)
%             barData(nUnit) = nanmean(nPeriodData(nUnit).unit_no_trial)*params.frameRate;
%         end
%         barSign     = -1;
%         barHistWithDist(barData(:), dist, '', barSeries, 'g', barSign); 
%         hold off;
%     end
% end
% 
% setPrint(8, 6, 'FR_distribution', 'pdf')