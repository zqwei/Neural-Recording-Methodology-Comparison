% 
% to first obtain the linear filtered signal
% to second obtain nonlinearity and noise

addpath('../Func');
setDir;
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
load([TempDatDir 'DataListS2CS1Model.mat']);

% nonlinear function
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));

int_noise = [3.0];
ext_noise = [0.45];

for nData = 1:1
    
    DataSetList(1+nData) = DataSetList(nData);
    DataSetList(1+nData).name = ['FineTuned' DataSetList(nData).name];
    
    load([TempDatDir DataSetList(nData).name '.mat'], 'nDataSet');
    for nUnit  = 1:length(nDataSet)
        param  = squeeze(nlParams(2, nUnit, :));
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
    save([TempDatDir DataSetList(1+nData).name '.mat'], 'nDataSet');     
end
 
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(1).params);       
sizeGroup = histcounts(unitGroup, 0:3);
disp(sizeGroup/sum(sizeGroup))
disp(sum(sizeGroup))

save([TempDatDir 'DataListS2CS1Model.mat'], 'DataSetList');