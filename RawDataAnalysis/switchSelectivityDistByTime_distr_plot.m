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

for nData      = [1 3 4]
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat']);
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat']);
    end
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
    sizeGroup = histcounts(unitGroup, 0:3);
    disp(sizeGroup)
    bootstat = bootstrp(1000,@(x)histcounts(x, 0:3)/length(x),unitGroup);
    disp(mean(bootstat))
    if nData == 1
        refEphys.mono = bootstat(:, 2);
        refEphys.multi = bootstat(:, 3);
    end
    if nData == 4
        performanceMat(1).mono = bootstat(:, 2);
        performanceMat(1).multi = bootstat(:, 3);
    end
    if nData == 3
        performanceMat(2).mono = bootstat(:, 2);
        performanceMat(2).multi = bootstat(:, 3);
    end
end
