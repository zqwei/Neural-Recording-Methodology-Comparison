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

load([TempDatDir DataSetList(3).name '.mat'])
nDataSet               = getFakeMLSpikeCaImagingData(nDataSet, params.frameRate);  
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 1;
DataSetList(nData).name    = 'ModelSpikeMLSpike_Ca_Slow_Short_Delay';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

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
load([TempDatDir DataSetList(4).name '.mat'])
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
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([TempDatDir 'DataListC2SRandomDeconvModel.mat'], 'DataSetList');

addpath('../Func');
setDir;
load ([TempDatDir 'DataListC2SRandomDeconvModel.mat']);
load([TempDatDir 'FineTuned6fNLParams.mat'], 'nlParams');
minNumTrialToAnalysis  = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
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

nlParams(sum(isnan(nlParams), 2)>0, :) = [];

load([TempDatDir 'Shuffle_Ca_Fast_SShort_Delay_withOLRemoval.mat'])
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
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
save([TempDatDir 'DataListC2SRandomDeconvModel.mat'], 'DataSetList');