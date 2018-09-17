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
load([TempDatDir 'FineTunedNLParamsFix.mat'], 'nlParams');
load([TempDatDir 'DataListS2CModelFix.mat']);

% nonlinear function
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));

int_noise = [2.0, 2.0, 0.95];
ext_noise = [0.15, 0.45, 0.40];

for nData = 1:3
    
    DataSetList(3+nData) = DataSetList(nData);
    DataSetList(3+nData).name = ['FineTuned' DataSetList(nData).name];
    param  = nanmedian(squeeze(nlParams(nData, :, :)));
    
    load([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
    for nUnit  = 1:length(nDataSet)
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
    save([TempDatDir DataSetList(3+nData).name '.mat'], 'nDataSet');     
end

save([TempDatDir 'DataListS2CModelFix.mat'], 'DataSetList');