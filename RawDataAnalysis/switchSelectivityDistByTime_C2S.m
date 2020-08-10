%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};
nDataList = [3 4 10];

figure;

for mData      = 1:3
    nData = nDataList(mData);
    load([TempDatDir DataSetList(nData).name '.mat']);
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    load([TempDatDir DataSetList(nData).name '_C2S MCMC model.mat']);
    unitGroup_ = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    for n = 1:3
        subplot(3, 6, n+mData*6-6)
        sizeGroup = histcounts(unitGroup_(unitGroup==n-1), 0:3);
        bar(sizeGroup/sum(sizeGroup))
        box off
        ylim([0, 1])
    end
    load([TempDatDir DataSetList(nData).name '_C2S MLSpike model.mat']);
    for n = 1:3
        subplot(3, 6, n+3+mData*6-6)
        sizeGroup = histcounts(unitGroup_(unitGroup==n-1), 0:3);
        bar(sizeGroup/sum(sizeGroup))
        box off
        ylim([0, 1])
    end
end


saveas(gcf, 'Selectivity_after_C2S.svg')