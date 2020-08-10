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
model_                 = '_S2C Linear model';
nData                  = nData + 1;
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP43_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP43_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
params.expression      = 'Transgentic';
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP43.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% short Ca slow virus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Slow_Short_Delay_Virus';
model_                 = '_S2C Linear model';
nData                  = nData + 1;
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV.mat'], 'nDataSet');
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

model_                 = '_S2C Linear model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP517_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP517_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP517.mat'], 'nDataSet');
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

model_                 = '_S2C Linear model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP43_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP43_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_ALM_Long_GP43.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GP43 short delay #10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_Fast_SShort_Delay';
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round((params.polein - 0.5) * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Transgentic';

model_                 = '_S2C Linear model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP517_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP517_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_GP517.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike short delay imaging S1 6s-AAV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_                  = 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus';
params.frameRate       =  29.68/2;
params.binsize         =  1/params.frameRate;
params.polein          =  -2.6;
params.poleout         =  -1.3;
minTimeToAnalysis      =  round(-3.1 * params.frameRate);
maxTimeToAnalysis      =  round(2.0 * params.frameRate);
params.timeWindowIndexRange  = minTimeToAnalysis : maxTimeToAnalysis;
params.timeSeries      = params.timeWindowIndexRange * params.binsize;
params.minNumTrialToAnalysis =  minNumTrialToAnalysis;
params.expression      = 'Virus';

model_                 = '_S2C Linear model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV_S1_S2C_Linear.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Hill function model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV_S1_S2C_Hill.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

model_                 = '_S2C Sigmoid model';
nData                  = nData + 1;
DataSetList(nData).name    = [name_ model_];
DataSetList(nData).params  = params;
load([tempDatOld 'FineTunedModeled_6s_AAV_S1.mat'], 'nDataSet');
nonActiveNeuronIndex   = findNonActiveNeurons(nDataSet, params);
DataSetList(nData).ActiveNeuronIndex = ~nonActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');


save([TempDatDir 'DataListS2CShuffle.mat'], 'DataSetList');
