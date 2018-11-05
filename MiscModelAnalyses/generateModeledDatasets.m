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
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
clear params;
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir 'Shuffle_Spikes_Nuo_Short_Delay.mat'], 'nDataSet');
params      = DataSetList(1).params;
spikeDataSet= nDataSet;
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
clear DataSetList

% Modeled 6s-AAV
params.Fm              = ones(length(nDataSet), 1) * 21.3240;
params.Ca0             = ones(length(nDataSet), 1) * median(nlParaMat(:, 1));
params.n               = ones(length(nDataSet), 1) * median(nlParaMat(:, 2));
median_r               = 0.0505;
median_d               = 1.7064;
params.tau_r           = ones(length(nDataSet), 1) * median_r;
params.tau_d           = ones(length(nDataSet), 1) * median_d;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 1;
DataSetList(nData).name    = 'Modeled_6s_AAV_fix';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

% Modeled GP4.3
% tau_r and tau_d are given randomly
median_r               = 0.0927;
median_d               = 1.2294;
params.tau_r           = ones(length(nDataSet), 1) * median_r;
params.tau_d           = ones(length(nDataSet), 1) * median_d;
params.intNoise        = 1.5;
params.extNoise        = 4.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 2;
DataSetList(nData).name    = 'Modeled_GP43_fix';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

depth       = [spikeDataSet.depth_in_um];
spikeDataSet= spikeDataSet(depth<471);
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
params.Fm              = ones(length(nDataSet), 1) * 21.3240;
params.Ca0             = ones(length(nDataSet), 1) * median(nlParaMat(:, 1));
params.n               = ones(length(nDataSet), 1) * median(nlParaMat(:, 2));
% Modeled GP5.17
median_r               = 0.0213;
median_d               = 0.5898;
params.tau_r           = ones(length(nDataSet), 1) * median_r;
params.tau_d           = ones(length(nDataSet), 1) * median_d;
params.intNoise        = 1.5;
params.extNoise        = 4.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 3;
DataSetList(nData).name    = 'Modeled_GP517_fix';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet'); 

save([TempDatDir 'DataListS2CModelFix.mat'], 'DataSetList');