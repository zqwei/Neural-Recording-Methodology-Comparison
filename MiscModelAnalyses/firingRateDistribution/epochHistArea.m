% Examine the firing rates in ephys vs whole cell recording data

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

% ephys

figure;
barSeries   = 0:0.1:5;
time1       = 3;
time2       = 4;

for nData = 1
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    timePoints  = [timePoints(time1), timePoints(time2)];
    load([TempDatDir DataSetList(nData).name '.mat'])
    m        = 2;
    
    dist        = 'Poisson';
    nPlot       = 1;
    nPeriodData_pre = dataInPeriods(nDataSet, [1 timePoints(1)], nPlot);
    barDataMean = nan(length(nPeriodData_pre), 1);
    for nUnit   = 1:length(nPeriodData_pre)
        barDataMean(nUnit) = mean([nPeriodData_pre(nUnit).unit_yes_trial; nPeriodData_pre(nUnit).unit_no_trial]);
    end
    
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
    hold on;
    barData     = nan(length(nPeriodData), 2);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit,1) = mean(nPeriodData(nUnit).unit_yes_trial);%max(mean(nPeriodData(nUnit).unit_yes_trial)/barDataMean(nUnit)-1, 0);
    end
    barSign     = 1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'b', barSign); 
%     barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 2) = mean(nPeriodData(nUnit).unit_no_trial);%max(mean(nPeriodData(nUnit).unit_no_trial)/barDataMean(nUnit)-1, 0);
    end
%     barSign     = -1;
    barData(barData==0) = nan;
    barHistWithDist(nanmean(barData, 2), dist, '', barSeries, 'r', barSign); 
    hold off;
end

disp(nanmean(nanmean(barData, 2)))
disp(nanstd(nanmean(barData, 2)))

% ca_factor = [9.5005, 1.9279, 5.4207];
% 
% ca_factor = ca_factor./[60, 55, 55] * 15;

ca_factor = [2.3701, 0.3649, 1.2353];

hold on;
for nData = 4
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    timePoints  = [timePoints(time1), timePoints(time2)];
    load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    m        = 2;
    
    dist        = 'Poisson';
    nPlot       = 1;
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
    hold on;
    barData     = nan(length(nPeriodData), 2);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 1) = sum(max(nPeriodData(nUnit).unit_yes_trial, 0))/ca_factor(1)/(params.timeSeries(timePoints(2))-params.timeSeries(timePoints(1)));
    end
%     barSign     = 1;
%     barHistWithDist(barData(:), dist, '', barSeries, 'g', barSign); 
%     barData     = nan(length(nPeriodData), 1);
    for nUnit   = 1:length(nPeriodData)
        barData(nUnit, 2) = sum(max(nPeriodData(nUnit).unit_no_trial, 0))/ca_factor(1)/(params.timeSeries(timePoints(2))-params.timeSeries(timePoints(1)));
    end
%     barSign     = -1;
    barData(barData==0) = nan;
    barHistWithDist(nanmean(barData, 2), dist, '', barSeries, 'k', barSign); 
    hold off;
end

disp(nanmean(nanmean(barData, 2)))
disp(nanstd(nanmean(barData, 2)))

xlabel('(R-R_{pre})/R_{pre}')
ylabel('Frequency')
legend({'Ephys', '6s-AAV'})
box off
set(gca, 'TickDir', 'out')


load ([TempDatDir 'DataListC2SMCMCSingleTrialModel.mat']);
nDatalist= 1;
barSeries   = 0:1:15;
for nData = nDatalist
    params   = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    timePoints  = [timePoints(time1), timePoints(time2)];
    load([TempDatDir DataSetList(nData).name '.mat'])
    m        = 2;

    dist        = 'Poisson';
    nPlot       = 1;
    nPeriodData = dataInPeriods(nDataSet, timePoints, nPlot);
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

disp(nanmean(nanmean(barData, 2)))
disp(nanstd(nanmean(barData, 2)))

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