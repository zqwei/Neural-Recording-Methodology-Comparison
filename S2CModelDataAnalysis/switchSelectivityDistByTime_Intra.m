%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CIntraModel.mat']);

cmap = [0.8000    0.8000    0.8000;
       1.0000    0.6000         0;
       0    0.8000         0];

for nData      = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])    
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);       
    sizeGroup = histcounts(unitGroup, 0:3);
    % disp([sizeGroup(2), sizeGroup(3)])
    figure('Visible', 'off');
    groupNames      = {'Non.', 'Homo.', 'Dynamical'};
    pie(sizeGroup)
    disp(sizeGroup/sum(sizeGroup))
    disp(sum(sizeGroup))
%     colormap(cmap)
%     set(gca, 'TickDir', 'out')
%     setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name])
end