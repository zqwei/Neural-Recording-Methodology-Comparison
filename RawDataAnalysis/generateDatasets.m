%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir 'ALM_compiled_all_data.mat'],'cellType_all')


minNumTrialToAnalysis  = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay Nuo #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5; % Hz per epoch
nDataSet               = getSpikeDataWithEphysTime(SpikingShortNuoDir, SpikingShortNuoFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize, cellType_all);
DataSetList(1).name    = 'Shuffle_Spikes_Nuo_Short_Delay';
DataSetList(1).params  = params; 
DataSetList(1).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(1).name '_old.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike -- DF/F (removed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nDataSet               = getDFFSpike(nDataSet, params);
% DataSetList(5).name    = 'Shuffle_DFFSpikes';
% DataSetList(5).params  = params; 
% DataSetList(5).ActiveNeuronIndex = ~findNonActiveNeurons(nDataSet, params) & DataSetList(1).ActiveNeuronIndex;
% save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca fast #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.4;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingData(CaImagingShortDelayFastDir, ...
                                          CaImagingShortDelayFastFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(2).name    = 'Shuffle_Ca_Fast_Short_Delay';
DataSetList(2).params  = params; 
DataSetList(2).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(2).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.4;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingData(CaImagingShortDelaySlowDir, ...
                                          CaImagingShortDelaySlowFileList, ...
                                          params.minNumTrialToAnalysis, params); 
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(3).name    = 'Shuffle_Ca_Slow_Short_Delay';
DataSetList(3).params  = params; 
DataSetList(3).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(3).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus #4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.4;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Virus';
nDataSet               = getCaImagingData(CaImagingShortDelaySlowVirusDir, ...
                                          CaImagingShortDelaySlowVirusFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(4).name    = 'Shuffle_Ca_Slow_Short_Delay_Virus';
DataSetList(4).params  = params; 
DataSetList(4).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(4).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca fast #5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.poleout         =  -3.0;
params.polein          =  params.poleout -1.2;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingData(CaImagingLongDelayFastDir, ...
                                          CaImagingLongDelayFastFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(5).name    = 'Shuffle_Ca_Fast_Long_Delay';
DataSetList(5).params  = params; 
DataSetList(5).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(5).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca slow #6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.poleout         =  -3.0;
params.polein          =  params.poleout -1.2;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingData(CaImagingLongDelaySlowDir, ...
                                          CaImagingLongDelaySlowFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(6).name    = 'Shuffle_Ca_Slow_Long_Delay';
DataSetList(6).params  = params; 
DataSetList(6).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(6).name '.mat'], 'nDataSet');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike long delay Nuo #7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.poleout         =  -3.0;
params.polein          =  params.poleout -1.2;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5; % Hz per epoch
nDataSet               = getSpikeDataWithEphysTime(SpikingLongNuoDir, SpikingLongNuoFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize, []);
DataSetList(7).name    = 'Shuffle_Spikes_Nuo_Long_Delay';
DataSetList(7).params  = params; 
DataSetList(7).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(7).name '.mat'], 'nDataSet'); % no cell depth info



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay Hi #8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5; % Hz per epoch
nDataSet               = getSpikeHiDataWithEphysTime(SpikingShortHiDir, SpikingShortHiFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize);
DataSetList(8).name    = 'Shuffle_Spikes_Hi_Short_Delay';
DataSetList(8).params  = params; 
DataSetList(8).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(8).name '_old.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay Hi intra #9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5; % Hz per epoch
nDataSet               = getHiIntraSpikeDataWithEphysTime(SpikingShortHiIntraDir, SpikingShortHiIntraFileList, params.timeSeries, params.binsize);
DataSetList(9).name    = 'Shuffle_Spikes_Hi_intra_Short_Delay';
DataSetList(9).params  = params; 
DataSetList(9).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(9).name '_old.mat'], 'nDataSet');


minNumTrialToAnalysis  = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP43 short delay #10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  30.0255/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.4;
params.poleout         =  -1.2;
start_time             =  params.polein - 0.5;
minTimeToAnalysis      =  round(start_time * params.frameRate);
maxTimeToAnalysis      =  round(1.2 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
nDataSet               = getCaImagingDataKD(CaImagingSShortDelayFastDir, ...
                                          CaImagingSShortDelayFastFileList, ...
                                          params.minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(10).name    = 'Shuffle_Ca_Fast_SShort_Delay';
DataSetList(10).params  = params; 
DataSetList(10).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(10).name '.mat'], 'nDataSet');


minNumTrialToAnalysis  = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay S1 #11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  7;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'None';
minFiringRate          = 5; % Hz per epoch
nDataSet               = getSpikeDataWithEphysTime(SpikingS1NuoDir, SpikingS1NuoFileList, params.minNumTrialToAnalysis, params.timeSeries, params.binsize, []);
DataSetList(11).name    = 'Shuffle_Spikes_Nuo_S1_Short_Delay';
DataSetList(11).params  = params; 
DataSetList(11).ActiveNeuronIndex = findHighFiringUnits(nDataSet, params, minFiringRate);
save([TempDatDir DataSetList(11).name '_old.mat'], 'nDataSet');


minNumTrialToAnalysis  = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay imaging S1 #12 6s-AAV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.frameRate       =  7;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.07;
params.poleout         =  -1.05;
minTimeToAnalysis      =  round(-2.6 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Virus';
nDataSet               = getCaImagingDataSP(CaImagingS1SlowVirusDir, CaImagingS1SlowVirusFileList, minNumTrialToAnalysis, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(12).name    = 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus';
DataSetList(12).params  = params; 
DataSetList(12).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(12).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileList            = {SpikingShortNuoFileList; ...
                       CaImagingShortDelayFastFileList; ...
                       CaImagingShortDelaySlowFileList;...
                       CaImagingShortDelaySlowVirusFileList;...
                       CaImagingLongDelayFastFileList;...
                       CaImagingLongDelaySlowFileList;...
                       SpikingLongNuoFileList;...
                       SpikingShortHiFileList;...
                       SpikingShortHiIntraFileList;...
                       CaImagingSShortDelayFastFileList;...
                       SpikingS1NuoFileList;...
                       CaImagingS1SlowVirusFileList};

for nData           = 1:length(fileList)
    if exist([TempDatDir DataSetList(nData).name '_old.mat'], 'file')
        load([TempDatDir DataSetList(nData).name '_old.mat']);
    else
        load([TempDatDir DataSetList(nData).name '.mat']);
    end
    DataSetList(nData).cellinfo  = repmat(struct('fileName',1, 'nUnit', 1, ...
                                'AP_axis',1, 'ML_axis', 1, 'depth', 1,...
                                'expression','none', 'cellType',-1),length(nDataSet), 1);
    DataSetList(nData).ROCIndex  = ROCPop(nDataSet, DataSetList(nData).params);
    for nUnit  = 1:length(nDataSet)
        nFileList                                     = fileList{nData};
        if nData ~= 12
            sessionIndex = nDataSet(nUnit).sessionIndex;
        else
            sessionIndex = floor(nDataSet(nUnit).sessionIndex/100);
        end
        DataSetList(nData).cellinfo(nUnit).fileName   = nFileList(sessionIndex).name;
        if isfield(nDataSet(nUnit), 'animalNameIndex')
            DataSetList(nData).cellinfo(nUnit).anmName= nDataSet(nUnit).animalNameIndex;
        else
            animalNameIndex                               = strfind(nFileList(sessionIndex).name, '_');
            DataSetList(nData).cellinfo(nUnit).anmName    = nFileList(sessionIndex).name(1:animalNameIndex(1));
        end
        DataSetList(nData).cellinfo(nUnit).nUnit      = nDataSet(nUnit).nUnit;
        DataSetList(nData).cellinfo(nUnit).AP_axis    = nDataSet(nUnit).AP_in_um;
        DataSetList(nData).cellinfo(nUnit).ML_axis    = nDataSet(nUnit).ML_in_um;
        DataSetList(nData).cellinfo(nUnit).depth      = nDataSet(nUnit).depth_in_um;
        DataSetList(nData).cellinfo(nUnit).expression = DataSetList(nData).params.expression;
        if strcmp(nDataSet(nUnit).cell_type, 'putative_interneuron')
            DataSetList(nData).cellinfo(nUnit).cellType   = 2;
        elseif strcmp(nDataSet(nUnit).cell_type, 'putative_pyramidal')
            DataSetList(nData).cellinfo(nUnit).cellType   = 1;
        elseif strcmp(nDataSet(nUnit).cell_type, 'unclassified')
            DataSetList(nData).cellinfo(nUnit).cellType   = 0;
        else
            DataSetList(nData).cellinfo(nUnit).cellType   = nDataSet(nUnit).cell_type;
        end        
    end
end


for nData     = [1 8:9]
    cellType = [DataSetList(nData).cellinfo.cellType];
    depth = [DataSetList(nData).cellinfo.depth];
    validCellIndex = cellType==1 & depth>100 & depth<800;
    DataSetList(nData).ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex(validCellIndex);
    DataSetList(nData).cellinfo = DataSetList(nData).cellinfo(validCellIndex);
    DataSetList(nData).ROCIndex = DataSetList(nData).ROCIndex(validCellIndex, :);
    load([TempDatDir DataSetList(nData).name '_old.mat'], 'nDataSet');
    nDataSet = nDataSet(validCellIndex);
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
end


for nData     = [11]
    cellType = [DataSetList(nData).cellinfo.cellType];
    depth = [DataSetList(nData).cellinfo.depth];
    validCellIndex = cellType==1;% & depth>100 & depth<450;
    DataSetList(nData).ActiveNeuronIndex = DataSetList(nData).ActiveNeuronIndex(validCellIndex);
    DataSetList(nData).cellinfo = DataSetList(nData).cellinfo(validCellIndex);
    DataSetList(nData).ROCIndex = DataSetList(nData).ROCIndex(validCellIndex, :);
    load([TempDatDir DataSetList(nData).name '_old.mat'], 'nDataSet');
    nDataSet = nDataSet(validCellIndex);
    save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
end

save([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
