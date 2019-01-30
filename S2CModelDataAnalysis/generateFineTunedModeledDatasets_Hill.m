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
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
load([TempDatDir 'FineTunedNLParams_Hill.mat'], 'nlParams');
DataSetList_name = {'Modeled_6s_AAV', ...
                     'Modeled_GP43', ...
                     'Modeled_GP517', ...
                     'Modeled_ALM_Long_GP43', ...
                     'Modeled_ALM_Long_GP517', ...
                     'Modeled_6s_AAV_S1'};

paramIndex = [1, 1, 1, 7, 7, 11];
% nonlinear function
g = @(p,x) p(1) + p(2)*x.^p(3)/(p(4)+x.^p(3));

int_noise = [4.0, 4.0, 4.0, 2.0, 3.0, 3.0];
ext_noise = [0, 0, 0, 0, 0, 0];

for nData = 1:length(DataSetList_name)
    DataSetListName = ['FineTuned' DataSetList_name{nData} '_S2C_Linear'];
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
    
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(7).params);       
    sizeGroup = histcounts(unitGroup, 0:3);
    disp(sizeGroup/sum(sizeGroup))
    disp(sum(sizeGroup))
    
    save([TempDatDir DataSetListName '.mat'], 'nDataSet');     
end