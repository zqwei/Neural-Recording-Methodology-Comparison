%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
% 
% Figure 2, Wei et al., 2020
%
% 
%
% author: Ziqiang Wei
% email: weiz@janelia.hhmi.org
%
% 
% 

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};

for nData      = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    sizeGroup = histcounts(unitGroup, 0:3);
    figure('Visible', 'off');
    groupNames      = {['Non' newline 'n = ' num2str(sum(unitGroup==0))], ... 
                       ['Mono' newline 'n = ' num2str(sum(unitGroup==1))], ...
                       ['Multi' newline 'n = ' num2str(sum(unitGroup==2))]};
    donut(sizeGroup, groupNames, groupColors);
    axis off
    legend('Location','eastoutside')
    legend('boxoff')
    setPrint(8, 6, [PlotDir 'WebsitePlots/' DataSetList(nData).name '_selectivity'], 'svg')
end

close all