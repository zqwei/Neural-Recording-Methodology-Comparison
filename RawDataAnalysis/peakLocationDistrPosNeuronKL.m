addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsPeakLocation'],'dir')
    mkdir([PlotDir 'SingleUnitsPeakLocation'])
end

for nData     = [1 3 4 10]%1:length(DataSetList)%[11 12]%
    if ~exist([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '.mat'])
        neuronRemoveList = false(length(nDataSet), 1);
    else
        load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
    
    disp(DataSetList(nData).name)
    params      = DataSetList(nData).params;
    timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
    
    numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
    yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
    noProfileMatrix     = yesProfileMatrix;
    positivePeak        = false(length(nDataSet));
    for nUnit        = 1:length(nDataSet)
        yesData      = mean(nDataSet(nUnit).unit_yes_trial);
        noData       = mean(nDataSet(nUnit).unit_no_trial);
        maxData      = max([yesData, noData]);
        minData      = min([yesData, noData]);
        rData        = (maxData - minData);
        yesData      = (yesData - minData)/(maxData - minData);
        noData       = (noData - minData)/(maxData - minData);
        yesProfileMatrix(nUnit, :)    = yesData;
        noProfileMatrix(nUnit, :)     = noData;
        positivePeak(nUnit)       = mean(yesData(timePoints(1):timePoints(2))) <= mean(yesData(timePoints(2):timePoints(4))) ...
                                   || mean(noData(timePoints(1):timePoints(2))) <= mean(noData(timePoints(2):timePoints(4)));

    end
    actMat        = [yesProfileMatrix, noProfileMatrix];
    actMat        = actMat(positivePeak, :);
    [~, maxId]    = max(actMat, [], 2);

    timeStep  = DataSetList(nData).params.timeSeries;
    timeTag   = timePoints(2):timePoints(4)+13; % sample to response
    numTime   = length(timeTag);
    polein    = DataSetList(nData).params.polein;
    poleout   = DataSetList(nData).params.poleout;
    kl        = zeros(nNeuron, 1);
    for nNeuron    = 1:length(maxId)
        maxId_     = maxId;
        maxId_(nNeuron) = [];
        countMaxId = (hist(maxId_, 1:numTimeBin*2)+1)/size(actMat,1);
        countMaxId = countMaxId / sum(countMaxId);
        % kl(nNeuron)= 1/length(countMaxId)*sum(log(countMaxId*length(countMaxId)));
        [emd_, emd_val] = emd_distance((1:length(countMaxId))', (1:length(countMaxId))', (ones(size(countMaxId))/length(countMaxId))', countMaxId'/sum(countMaxId), @gdf);
        kl(nNeuron)= emd_val;
    end
    disp([mean(kl) std(kl)])
end

close all



% % % % Shuffle_Spikes_Nuo_Short_Delay
% % % %     5.5408    0.0183
% % % % 
% % % % Shuffle_Ca_Slow_Short_Delay
% % % %     6.2397    0.0198
% % % % 
% % % % Shuffle_Ca_Slow_Short_Delay_Virus
% % % %     3.0310    2.5735
% % % % 
% % % % Shuffle_Ca_Fast_SShort_Delay
% % % %     5.1598    0.0101