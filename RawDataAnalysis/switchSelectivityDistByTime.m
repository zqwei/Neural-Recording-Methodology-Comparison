%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};

for nData      = [1 3 4 9 11 12]%1:length(DataSetList)
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    sizeGroup = histcounts(unitGroup, 0:3);
    disp(sizeGroup/sum(sizeGroup))
    disp(sum(sizeGroup))
    figure('Visible', 'off');
    groupNames      = {['Non' newline 'n = ' num2str(sum(unitGroup==0))], ... 
                       ['Mono' newline 'n = ' num2str(sum(unitGroup==1))], ...
                       ['Multi' newline 'n = ' num2str(sum(unitGroup==2))]};
    donut(sizeGroup, groupNames, groupColors);
    % pie(sizeGroup, groupNames)
    axis off
    legend('Location','eastoutside')
    legend('boxoff')
    setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTime_' DataSetList(nData).name], 'png')
%     setPrint(8, 6, [PlotDir 'SingleUnitsTscore/' DataSetList(nData).name '_selectivity'], 'svg')
% 
%     depth                        = [DataSetList(nData).cellinfo(:).depth];
%     depth                        = depth(~neuronRemoveList);
%     depthStart                   = 100;
%     depthBin                     = 50;
%     depthEnd                     = 900;
%     depth                        = floor((depth-depthStart)/depthBin)*depthBin+depthStart;
%     depth(depth>depthEnd)        = depthEnd;
%     depth(depth<depthStart)      = depthStart;
% 
%     uniqueDepth                  = depthStart:depthBin:depthEnd;
%     depthStrings                 = cell(length(uniqueDepth),1);
%     depthStrings(1:end)          = {''};
%     if length(uniqueDepth)       <=3
%         depthStrings             = cellstr(num2str(uniqueDepth'));
%     else
%         stepLength               = floor(length(uniqueDepth)/3);
%         depthStrings(1:stepLength:end) = cellstr(num2str(uniqueDepth(1:stepLength:end)'));
%     end
% 
%     groupCounts = zeros(length(uniqueDepth), 3);
% 
%     for nGroup = 1:3
%         nUnitGroup = unitGroup == nGroup-1;
%         for nDepth = 1:length(uniqueDepth)
%             groupCounts(nDepth, nGroup) = sum(nUnitGroup & depth' == uniqueDepth(nDepth));
%         end
%     end
% 
%     [tab, chi2, p] = crosstab(unitGroup(unitGroup>0), depth(unitGroup>0)');
%     disp(p)
% 
%     groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
% 
%     figure('Visible', 'off');
%     subplot(1, 2, 1)
%     barh(-uniqueDepth, groupPerCounts, 'stack', 'edgecolor', 'none');
%     xlim([0 1])
%     box off
%     xlabel('% cell type')
%     ylabel('Depth (um)')
%     set(gca, 'TickDir', 'out')
%     ylim([-850 0])
%     set(gca, 'yTick', -800:400:0)
%     colormap(cmap)
% 
%     subplot(1, 2, 2)
%     barh(-uniqueDepth, sum(groupCounts,2),'k')
%     xlabel('# cells')
%     ylabel('Depth (um)')
%     ylim([-850 0])
%     set(gca, 'yTick', -800:400:0)
%     set(gca, 'TickDir', 'out')
%     setPrint(8*2, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreTimeDepth_' DataSetList(nData).name])
end

% % % % figure;
% % % % hold on
% % % % for nColor = 1:length(groupNames)
% % % %     plot(0, nColor, 's', 'color', cmap(nColor,:), 'MarkerFaceColor',cmap(nColor,:),'MarkerSize', 8)
% % % %     text(1, nColor, groupNames{nColor})
% % % % end
% % % % xlim([0 10])
% % % % hold off
% % % % axis off
% % % % setPrint(3, 2, [PlotDir 'SingleUnitsTscore/SingleUnitsTscore_Label'])
% % % % close all
% % % % 
% % % % for nData = [1 4]
% % % %     if nData   == 1
% % % %         load([TempDatDir DataSetList(nData).name '.mat'])
% % % %         neuronRemoveList = false(length(nDataSet), 1);
% % % %     else
% % % %         load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
% % % %     end
% % % %     [~, ~, anmIndex] = unique(cell2mat({DataSetList(nData).cellinfo.anmName}'), 'rows');
% % % %     numGroup    = 3;
% % % %     groupCounts = zeros(anmIndex(end), numGroup);
% % % %     anmIndex  = anmIndex(~neuronRemoveList);
% % % %     unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
% % % %     for nAnm    = 1:anmIndex(end)
% % % %         if sum(anmIndex == nAnm) > 50
% % % %             for nGroup = 1:numGroup
% % % %                 groupCounts(nAnm, nGroup) = sum(unitGroup == nGroup-1 & anmIndex == nAnm);
% % % %             end
% % % %         end
% % % %     end
% % % %     zeroGroups  = sum(groupCounts, 2) == 0;
% % % %     groupCounts = groupCounts(~zeroGroups, :);
% % % % %     [tab, chi2, p] = crosstab(unitGroup(unitGroup>0), anmIndex(unitGroup>0)');
% % % %     groupPerCounts = bsxfun(@rdivide, groupCounts, sum(groupCounts, 2));
% % % %     figure
% % % %     subplot(2, 1, 1)
% % % %     barh(groupPerCounts, 'stack', 'edgecolor', 'none');
% % % %     colormap(cmap)
% % % %     %     caxis([1 8])
% % % %     xlim([0 1])
% % % %     box off
% % % %     xlabel('% cell type')
% % % %     ylabel('Animal index')
% % % %     set(gca, 'yTickLabel', {})
% % % %     set(gca, 'TickDir', 'out')
% % % %     subplot(2, 1, 2)
% % % %     barh(sum(groupCounts,2), 'k')
% % % %     xlabel('# cells')
% % % %     ylabel('Animal index')
% % % %     set(gca, 'yTickLabel', {})
% % % %     set(gca, 'TickDir', 'out')
% % % %     setPrint(8, 6*2, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreANM_' DataSetList(nData).name])
% % % % end
% % % % 
% % % % load ([TempDatDir 'DataListShuffleConfounding.mat']);
% % % % nData          = 5;
% % % % load([TempDatDir DataSetList(nData).name '.mat'])
% % % % unitGroup      = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
% % % % figure;
% % % % areaIndex      = false(length(unitGroup), 1);
% % % % sizeGroup      = histcounts(unitGroup, 0:3);
% % % % pie(sizeGroup)
% % % % colormap(cmap)
% % % % disp(sum(sizeGroup))
% % % % setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreAPML_Pop_' DataSetList(nData).name])
% % % % figure;
% % % % APLoc          = [DataSetList(nData).cellinfo.AP_axis];
% % % % MLLoc          = [DataSetList(nData).cellinfo.ML_axis];
% % % % sizeGroup      = histcounts(unitGroup(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900), 0:3);
% % % % areaIndex(APLoc>2400 & APLoc<2600 & MLLoc>1100 & MLLoc<1900) = true;
% % % % [tab, chi2, p] = crosstab(unitGroup(unitGroup>0), areaIndex(unitGroup>0)');
% % % % groupNames     = {'Non.', 'Homo.', 'Dynamicial'};
% % % % pie(sizeGroup)
% % % % colormap(cmap)
% % % % disp(sum(sizeGroup))
% % % % setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreAPML_Sub_' DataSetList(nData).name])
% % % % 
% % % % figure
% % % % hold on
% % % % plot(MLLoc, APLoc, '.', 'color', [0.5 0.5 0.5])
% % % % plot(MLLoc(areaIndex), APLoc(areaIndex), '.k')
% % % % setPrint(8, 6, [PlotDir 'SingleUnitsTscore/SingleUnitsTscoreAPML_Sub_Loc_' DataSetList(nData).name])
