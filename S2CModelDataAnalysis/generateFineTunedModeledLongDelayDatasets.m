%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating new data using fine tuned params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% one should run:
% 1.
% generateModeledDatasets
% 2.
% nonlinearityTuning
% 
% to first obtain the linear filtered signal
% to second obtain nonlinearity and noise

addpath('../Func');
setDir;
load([TempDatDir 'FineTunedNLParamsLongDelay.mat'], 'nlParams');
DataSetList_name = {'Modeled_ALM_Long_GP43', 'Modeled_ALM_Long_GP517'};

% nonlinear function
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));

int_noise = [0.3, 0.5];
ext_noise = [0.15, 0.45];

for nData = 1:2
    DataSetListName = ['FineTuned' DataSetList_name{nData}];
    load([TempDatDir DataSetList_name{nData} '.mat'], 'nDataSet');
    for nUnit  = 1:length(nDataSet)
        param  = squeeze(nlParams(nData, nUnit, :));
        yesNoise = randn(size(nDataSet(nUnit).unit_yes_trial))*ext_noise(nData);
        noNoise  = randn(size(nDataSet(nUnit).unit_no_trial))*ext_noise(nData);
        yesIntNoise = randn(size(nDataSet(nUnit).unit_yes_trial))*int_noise(nData);
        noIntNoise  = randn(size(nDataSet(nUnit).unit_no_trial))*int_noise(nData);
        unit_yes_trial_linear = nDataSet(nUnit).unit_yes_trial_linear + yesIntNoise;
        unit_no_trial_linear  = nDataSet(nUnit).unit_no_trial_linear + noIntNoise;
        unit_yes_trial_linear(unit_yes_trial_linear<0) = 0;
        unit_no_trial_linear(unit_no_trial_linear<0)   = 0;
        
        nDataSet(nUnit).unit_yes_trial = g(param, unit_yes_trial_linear) + yesIntNoise;
        nDataSet(nUnit).unit_no_trial  = g(param, unit_no_trial_linear) + noIntNoise;
    end
    save([TempDatDir DataSetListName '.mat'], 'nDataSet');     
end

load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');

for nData      = 1:2
    DataSetListName = ['FineTuned' DataSetList_name{nData}];
    load([TempDatDir DataSetListName '.mat'])    
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(7).params);       
    sizeGroup = histcounts(unitGroup, 0:3);
    disp(sizeGroup/sum(sizeGroup))
    disp(sum(sizeGroup))
end