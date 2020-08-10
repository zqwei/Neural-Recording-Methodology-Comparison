%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Generate raw activity of all neurons sorted in different ways
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_V2_0

addpath('../Func');
setDir;
tempDatOld = '../TempDat_2019_01_28/';
nData = 0;

minNumTrialToAnalysis  = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Slow_Short_Delay';
model_                 = '_C2S MCMC model';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_Slow_Short_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Slow_Short_Delay_Virus';
model_                 = '_C2S MCMC model';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
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
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long Ca fast #5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Fast_Long_Delay';
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

model_                 = '_C2S MCMC model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_Fast_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_Fast_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % long Ca slow #6
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Slow_Long_Delay';
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.poleout         =  -3.0;
params.polein          =  params.poleout -1.2;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;

model_                 = '_C2S MCMC model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;

load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_Slow_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_Slow_Long_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP43 short delay #10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Fast_SShort_Delay';
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

model_                 = '_C2S MCMC model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_Fast_SShort_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_Fast_SShort_Delay.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay imaging S1 6s-AAV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus';
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

model_                 = '_C2S MCMC model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_SP_S1_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S MLSpike model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeMLSpike_Deconv_Ca_SP_S1_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_C2S Random linear deconv';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'ModelSpikeRandom_Deconv_Ca_SP_S1_Slow_Short_Delay_Virus.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


save([TempDatDir 'DataListC2SShuffle.mat'], 'DataSetList');
