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
load([TempDatDir 'Shuffle_Spikes_Nuo_S1_Short_Delay.mat'], 'nDataSet');
params      = DataSetList(1).params;
spikeDataSet           = nDataSet;
ActiveNeuronIndex = DataSetList(1).ActiveNeuronIndex;
clear DataSetList

% Modeled 6s-AAV
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
randPairs              = randi([1 length(nlParaMat)], length(nDataSet), 1);
params.Ca0             = nlParaMat(randPairs, 1);
params.n               = nlParaMat(randPairs, 2);
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingDataSigmoid(spikeDataSet, params);
nData                      = 1;
DataSetList(nData).name    = 'Modeled_6s_AAV_S1';
DataSetList(nData).params  = params; 
DataSetList(nData).ActiveNeuronIndex = ActiveNeuronIndex;
save([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');

save([TempDatDir 'DataListS2CS1Model.mat'], 'DataSetList');