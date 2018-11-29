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

minNumTrialToAnalysis  = 15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
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


load([TempDatDir 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus_withOLRemoval.mat'])

nDataSet               = getFakeSpikeSingleTrialOOPSIData(nDataSet);  
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
nData                  = 1;
DataSetList(nData).name    = 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_SP_S1_Slow_Short_Delay_Virus';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

save([TempDatDir 'DataListC2SMCMCSingleTrialS1Model.mat'], 'DataSetList');
