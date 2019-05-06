addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsContraIpsi'],'dir')
    mkdir([PlotDir 'SingleUnitsContraIpsi'])
end

cmap = cbrewer('qual', 'Set1', 3, 'cubic');

nIndex = 0;
figure;
hold on
for nData = [4 1 13]%[5 2 14] %[1 4 13]
    load([TempDatDir DataSetList(nData).name '.mat'])
    logPValueEpoch= getLogPValueTscoreSpikeEpoch(nDataSet, DataSetList(nData).params);
    unitGroup = plotTtestLogPSpikeEpoch (logPValueEpoch);
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
        contraIndex(nUnit)   = sum(noTrial(timePoints(2):end))<sum(yesTrial(timePoints(2):end));
    end
    contraCount = sum(unitGroup~=0 & contraIndex);
    ipsiCount = sum(unitGroup~=0 & ~contraIndex);
    nIndex = nIndex+1;
    bar(nIndex, contraCount/(contraCount+ipsiCount), ...
        'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
end

set(gca, 'xTick', 1:3)
set(gca, 'TickDir', 'out')
ylim([0.3 0.7])
xlim([0.5 nIndex+0.5])
setPrint(4, 5, 'Faction_mcmc', 'pdf')