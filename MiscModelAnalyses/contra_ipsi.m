addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat']);
load('temp_results.mat')
load('temp_results_var.mat')

figure;
nDataInd = [10 3 4];

for mData      = 1:3
    nData = nDataInd(mData);
    load([TempDatDir DataSetList(nData).name '.mat'])
    
    snr = (data_var{mData}.^2+1)./(data_snr{mData}.^2+1);
    p = [1 5:5:95 99];
    snr_thres = prctile(snr, p);
    
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params); 
    params      = DataSetList(nData).params;
    contraIndex = false(length(nDataSet), 1);
    yesActMat   = nan(length(nDataSet), length(params.timeSeries));
    noActMat    = nan(length(nDataSet), length(params.timeSeries));
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

    for nUnit   = 1:length(nDataSet)
        yesTrial = mean(nDataSet(nUnit).unit_yes_trial);
        noTrial  = mean(nDataSet(nUnit).unit_no_trial);
        yesActMat(nUnit, :)  = yesTrial;
        noActMat(nUnit, :)   = noTrial;
        if nData ==1 
            contraIndex(nUnit)   = sum(noTrial(timePoints(2):timePoints(5)))<sum(yesTrial(timePoints(2):timePoints(5)));
        else
            contraIndex(nUnit)   = sum(noTrial(timePoints(2):timePoints(4)))<sum(yesTrial(timePoints(2):timePoints(4)));    
        end
    end
    
    contra_ = zeros(20, 1);
    for n = 1:length(p)
        contraCount = (unitGroup~=0) & contraIndex & (snr>snr_thres(n));
        ipsiCount = (unitGroup~=0) & ~contraIndex & (snr>snr_thres(n));
        [sum(contraCount), sum(ipsiCount)]
%         contra_(n) = sum(contraCount)/(sum(contraCount)+sum(ipsiCount));
    end
%     
%     subplot(1, 3, mData)
%     hold on
%     plot(p, contra_, '-o')
%     xlim([0, 100])
%     ylim([0.3, 0.7])
%     xlabel('Percentile of SNR(%)')
%     ylabel('Fraction of cells')

    subplot(1, 3, mData)
    hold on
    [sorted_snr, p_snr] = sort(snr);
    n_ind = zeros(length(p_snr), 1);
    for n = 1:length(p_snr)
        n_ind(p_snr(n)) = n;
    end
    n_ind = n_ind/max(n_ind)*100;
    histogram(n_ind(unitGroup~=0 & contraIndex), 0:5:100)
    histogram(n_ind(unitGroup~=0 & ~contraIndex), 0:5:100)
    
end
% 
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperPosition', [0 0 5 1.5]);
% saveas(gcf, 'contra_v2.svg')