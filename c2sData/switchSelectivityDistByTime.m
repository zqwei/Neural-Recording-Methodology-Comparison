%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SShuffle.mat']);

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};

ffactor = [  0     0     1.5
             0     0     1.2
             0     0     1.5
             0     0     1.5
             0     0     1.5
             0     0     1.0]';

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    if contains(DataSetList(nData).name, 'Random')
        unitGroup = getLogPValueTscoreSpikeTimeAve(nDataSet, DataSetList(nData).params, ffactor(nData));
    else
        unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    end
    sizeGroup = histcounts(unitGroup, 0:3);
    figure('Visible', 'off');
    groupNames      = {['Non' newline 'n = ' num2str(sum(unitGroup==0))], ...
                       ['Mono' newline 'n = ' num2str(sum(unitGroup==1))], ...
                       ['Multi' newline 'n = ' num2str(sum(unitGroup==2))]};
    donut(sizeGroup, groupNames, groupColors);
    if contains(DataSetList(nData).name, 'Random')
        disp(DataSetList(nData).name)
        disp(sizeGroup/sum(sizeGroup))
    end
    axis off
    legend('Location','eastoutside')
    legend('boxoff')
    setPrint(8, 6, [PlotDir 'WebsitePlots/' DataSetList(nData).name '_selectivity'], 'svg')
end

close all
