%
% Compute distribution of selectivity for each single neuron
% 
% Here we drop the test of parameter of tau_r
% 
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

% all data is precomputed by code meanSynthEphysTraces.m

indexDatasets = [3, 4, 10];
% 2: short Ca GP517
% 3: short Ca slow
% 4: short Ca slow virus
% 5: long Ca fast
% 6: long Ca slow

for nData     = indexDatasets    
    params           = DataSetList(nData).params;
    sigma            = 0.15 / params.binsize; % 300 ms
    filterLength     = 11;
    filterStep       = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse      = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse      = filterInUse / sum (filterInUse); 

    per_list           = 0.02:0.01:0.98;
    noise_factor_list  = sqrt(icdf('Exponential', per_list, 0.3527));
    timeTag            = 8:60;
    
    numT      = length(params.timeSeries);
    yesData   = nan(97, length(DataSetList(nData).ActiveNeuronIndex), numT);
    noData    = nan(97, length(DataSetList(nData).ActiveNeuronIndex), numT);
    
    for nTau  = 1:97
        if exist([TempDatDir 'directNLDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'], 'file')
            load([TempDatDir 'directNLDeconv/' DataSetList(nData).name '_Tau' num2str(nTau, '%02d') '.mat'],'spikeDataSet')
            for nUnit        = 1:length(spikeDataSet)
                nUnitData    = spikeDataSet(nUnit).unit_yes_trial;
                yesUnitData  = getGaussianPSTH (filterInUse, nUnitData, 2);
                nUnitData    = spikeDataSet(nUnit).unit_no_trial;
                noUnitData   = getGaussianPSTH (filterInUse, nUnitData, 2);
                mean_rate    = min([yesUnitData, noUnitData]);
                yesUnitData  = yesUnitData - mean_rate;
                noUnitData   = noUnitData - mean_rate;
                yesData(nTau, nUnit, : ) = yesUnitData;
                noData(nTau, nUnit, :)   = noUnitData;
            end
            clear spikeDataSet
        end
    end
    
    save([TempDatDir 'directNLDeconv_' DataSetList(nData).name '.mat'], 'yesData', 'noData', '-V7.3');
    
end