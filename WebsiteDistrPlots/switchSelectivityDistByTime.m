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
x_labels = {'Fraction of mono-cells', 'Fraction of multi-cells'};
xi = 0.0:0.01:1.0; 

fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);
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
            set(gca,'XTick',[0, 1], 'YTick', [0, 1], 'TickDir', 'out')
            ylabel('Prob. density')
            xlabel(x_labels{n_})
            box off
            set(gca,'fontsize', 14)
            setPrint(8,4.5,[Result_ DataSetList(nData).name '_' file_ext{n_}], 'svg')
            close all;
        end
    end
end