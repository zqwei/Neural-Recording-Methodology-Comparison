%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
% load ([TempDatDir 'DataListShuffle.mat']);
% load ([TempDatDir 'DataListC2SShuffle.mat']);
load ([TempDatDir 'DataListS2CShuffle.mat']);
file_ext = {'MonoCell', 'MultiCell'};
xi = 0.0:0.01:1.0; 

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    sizeGroup = histcounts(unitGroup, 0:3);
    bootstat = bootstrp(1000,@(x)histcounts(x, 0:3)/length(x),unitGroup);
    for n_ = 1:2
        figure
        tmp_s = bootstat(:, n_+1);
        fs = ksdensity(tmp_s, xi);
        f_ = max(fs);
        plot(xi, fs/f_, '-k', 'linewid', 2);
        set(gca,'XTick',[0, 1], 'YTick', [], 'TickDir', 'out','Ycolor','none')
        box off
        setPrint(8,3,['../Plot/WebsiteDistrPlots/' DataSetList(nData).name '_' file_ext{n_}], 'svg')
        close all;
    end
end
