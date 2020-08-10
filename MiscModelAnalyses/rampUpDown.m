addpath('../Func');
setDir;
load('temp_results.mat')
load('temp_results_var.mat')
load([TempDatDir 'DataListShuffle.mat']);

figure;
nDataInd = [10 3 4];
for mData      = 1:3
    nData = nDataInd(mData);
    load([TempDatDir DataSetList(nData).name '.mat'])
    depth       = [nDataSet.depth_in_um];
    spikeDataSet= nDataSet;  
    snr = (data_var{mData}.^2+1)./(data_snr{mData}.^2+1);
    p = [1 5:5:95 99];
    snr_thres = prctile(snr, p);


    yesRampDown          = zeros(length(spikeDataSet), 1);
    noRampDown           = zeros(length(spikeDataSet), 1);
    timePoints(1)        = sum(DataSetList(nData).params.timeSeries<DataSetList(nData).params.polein);
    timePoints(2)        = sum(DataSetList(nData).params.timeSeries<0);
    contraIndex          = false(length(spikeDataSet), 1);
    unitGroup            = getLogPValueTscoreSpikeTime(spikeDataSet, DataSetList(nData).params); 

    for nUnit            = 1:length(spikeDataSet)
        meanPreSample    = [mean(spikeDataSet(nUnit).unit_yes_trial(:, 1:timePoints(1)), 2); mean(spikeDataSet(nUnit).unit_no_trial(:, 1:timePoints(1)), 2)];
        meanYesSample    = mean(spikeDataSet(nUnit).unit_yes_trial(:, timePoints(1):timePoints(2)), 2);
        meanNoSample     = mean(spikeDataSet(nUnit).unit_no_trial(:, timePoints(1):timePoints(2)), 2);
        if ttest2(meanPreSample, meanYesSample, 'tail', 'right')
            yesRampDown(nUnit) = 1;
        elseif ttest2(meanPreSample, meanYesSample, 'tail', 'left')
            yesRampDown(nUnit) = -1;
        end    
        if ttest2(meanPreSample, meanNoSample, 'tail', 'right')
            noRampDown(nUnit)  = 1;
        elseif ttest2(meanPreSample, meanNoSample, 'tail', 'left')
            noRampDown(nUnit)  = -1;
        end    
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        contraIndex(nUnit)= sum(noTrial(timePoints(1):end))<sum(yesTrial(timePoints(1):end));
    end

    RampDown = (yesRampDown==1 & noRampDown>=0) | (yesRampDown>=0 & noRampDown==1);% & unitGroup>0;
    RampUp   = (yesRampDown==-1 & noRampDown<=0) | (yesRampDown<=0 & noRampDown==-1);%  & unitGroup>0;
    
    up_ = zeros(20, 1);
    down_ = zeros(20, 1);
    for n = 1:length(p)
        down_(n) = mean(RampDown(snr>snr_thres(n)));
        up_(n) = mean(RampUp(snr>snr_thres(n)));
    end
    
    subplot(1, 3, mData)
    hold on
    plot(p, up_, '-o')
    plot(p, down_, '-o')
    xlim([0, 100])
    ylim([0, 0.6])
    xlabel('Percentile of SNR(%)')
    ylabel('Fraction of cells')
end

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 1.5]);
saveas(gcf, 'ramp_up_down_v2.svg')