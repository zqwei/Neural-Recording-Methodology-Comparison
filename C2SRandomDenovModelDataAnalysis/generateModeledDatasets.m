%
% Comparison based on single unit acitivity
% Generating Ca++ imaging data from ephys data using Tsai-Wen's model
%
% -------------------------------------------------------------------------
% version 1.0
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
TempDatDirOld = '../TempDat_2019_01_28/';

minNumTrialToAnalysis  = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
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

load([TempDatDirOld DataSetList(3).name '_withOLRemoval.mat'])
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
std_r                  = 0.0375; % 0; %0.0375;
median_r               = 0.0927;
std_d                  = 0.5374; % 0; %0.5374;
median_d               = 1.2294;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 1;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_Slow_Short_Delay';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus
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
load([TempDatDirOld DataSetList(4).name '_withOLRemoval.mat'])
std_r                  = 0.0246; % 0; %0.0246;
median_r               = 0.0505;
std_d                  = 0.4588; % 0; %0.4588;
median_d               = 1.7064;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 2;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_Slow_Short_Delay_Virus';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');

minNumTrialToAnalysis  = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca fast
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

load([TempDatDirOld 'Shuffle_Ca_Fast_SShort_Delay_withOLRemoval.mat'])
truncatedNormal        = truncate(makedist('Normal'), -0.9, 1.5);
std_r                  = 0.0222;
median_r               = 0.0213;
std_d                  = 0.5390;
median_d               = 0.5898;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 3;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_Fast_SShort_Delay';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short 6s-AAV S1
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
load([TempDatDirOld 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus_withOLRemoval.mat'])
std_r                  = 0.0246; % 0; %0.0246;
median_r               = 0.0505;
std_d                  = 0.4588; % 0; %0.4588;
median_d               = 1.7064;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 4;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_SP_S1_Slow_Short_Delay_Virus';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca fast
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
load([TempDatDirOld 'Shuffle_Ca_Fast_Long_Delay_withOLRemoval.mat'])
std_r                  = 0.0222;
median_r               = 0.0213;
std_d                  = 0.5390;
median_d               = 0.5898;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 5;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_Fast_Long_Delay';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.binsize         =  1/params.frameRate;
params.poleout         =  -3.0;
params.polein          =  params.poleout -1.2;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
load([TempDatDirOld 'Shuffle_Ca_Slow_Long_Delay_withOLRemoval.mat'])
std_r                  = 0.0375; % 0; %0.0375;
median_r               = 0.0927;
std_d                  = 0.5374; % 0; %0.5374;
median_d               = 1.2294;
tau_r                  = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
tau_d                  = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
nDataSet               = getFakeSpikeDeconvData(nDataSet, tau_r, tau_d, params);
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 6;
DataSetList(nData).name    = 'ModelSpikeRandom_Deconv_Ca_Slow_Long_Delay';
DataSetList(nData).params  = params;
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDirOld DataSetList(nData).name '.mat'], 'nDataSet');


save([TempDatDirOld 'DataListC2SRandomDeconvModel.mat'], 'DataSetList');
