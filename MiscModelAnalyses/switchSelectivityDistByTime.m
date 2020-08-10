%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
TempDatDir = '../TempDat/';
load ([TempDatDir 'DataListShuffle.mat']);
load('temp_results.mat')
load('temp_results_var.mat')

figure;
nDataInd = [10 3 4];
for mData      = 1:3
    nData = nDataInd(mData);
    snr = (data_var{mData}.^2+1)./(data_snr{mData}.^2+1);
    p = [1 5:5:95 99];
    snr_thres = prctile(snr, p);
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    mono_ = zeros(20, 1);
    multi_ = zeros(20, 1);
    for n = 1:length(p)
        unitGroup_ = unitGroup(snr>snr_thres(n));
        sizeGroup = histcounts(unitGroup_, 0:3);
        mono_(n) = sizeGroup(2)/sum(sizeGroup);
        multi_(n) = sizeGroup(3)/sum(sizeGroup);
    end
    
    subplot(1, 3, mData)
    hold on
    plot(p, mono_, '-o')
    plot(p, multi_, '-o')
    plot(p, 1-multi_-mono_, '-o')
    xlim([0, 100])
    ylim([0, 1])
    xlabel('Percentile of SNR(%)')
    ylabel('Fraction of cells')
end

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 1.5]);
saveas(gcf, 'selectivity.svg')