addpath('../Func');
setDir;
TempDatDir = '../Backups/TempDat_2019_01_28/';

numNeuron  = 50;
numTrial   = 30;
frameRate  = 14.84;
numTimeBin = 61;
timeSeries = (0:60)/frameRate;
peakTime   = zeros(1, numTimeBin);
peakTime(11) = 1;
ref_peak_list_model = nan(numNeuron, 4);
intNoise        = 0.1;
extNoise        = 0.0;
fix_spk_num     = 5;

%% model 1
model_names{1} = 'Ephys';
ephys      = sqrt(mean((peakTime - 1/numTimeBin).^2))*100;
ref_peak_list_model(:, 1) = ephys;

%% model 2
model_names{2} = 'Fix Calcium';
ref_spikeTimes = cell(numNeuron, numTrial);
for nTrial   = 1:numTrial
    spkTime  = sort(rand(fix_spk_num, 1))/frameRate + timeSeries(10);
    for nNeuron = 1:numNeuron
        ref_spikeTimes{nNeuron, nTrial} = spkTime;
    end
end

load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
Fm              = 21.3240;
Ca0             = median(nlParaMat(:, 1));
n               = median(nlParaMat(:, 2));
median_r        = 0.0505;
median_d        = 1.7064;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData          = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 2) = ref_peak_list;

%% model 3 -- spike time varies
model_names{3} = 'Spike times';
ref_spikeTimes_model1 = cell(numNeuron, numTrial);
for nTrial   = 1:numTrial
    for nNeuron = 1:numNeuron
        spkTime  = sort(rand(fix_spk_num, 1))/frameRate + timeSeries(10);
        ref_spikeTimes_model1{nNeuron, nTrial} = spkTime;
    end
end

load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
Fm              = 21.3240;
Ca0             = median(nlParaMat(:, 1));
n               = median(nlParaMat(:, 2));
median_r        = 0.0505;
median_d        = 1.7064;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model1(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 3) = ref_peak_list;

%% model 4 -- firing rate varies
model_names{4} = 'Firing rates';
ref_spikeTimes_model2 = cell(numNeuron, numTrial);
rate_        = ceil(fix_spk_num*frameRate)+(-25:25);
for nTrial   = 1:numTrial
    for nNeuron = 1:numNeuron
        r_   = rate_(nNeuron);
        spkTime  = cumsum(exprnd(1/r_, 30, 1));
        spkTime  = spkTime(spkTime<=1/frameRate)+ timeSeries(10);
        ref_spikeTimes_model2{nNeuron, nTrial} = spkTime;
    end
end
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
Fm              = 21.3240;
Ca0             = median(nlParaMat(:, 1));
n               = median(nlParaMat(:, 2));
median_r        = 0.0505;
median_d        = 1.7064;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model2(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 4) = ref_peak_list;

%% model 5 -- decay time
model_names{5} = 'Decay time';
ref_spikeTimes_model3 = ref_spikeTimes;
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
randPairs       = randi([1 length(nlParaMat)], numNeuron, 1);
truncatedNormal = truncate(makedist('Normal'), -1.5, 1.5);
Fm              = 21.3240;
Ca0             = median(nlParaMat(:, 1));
n               = median(nlParaMat(:, 2));
std_d           = 0.4588;
median_d        = 1.7064;
tau_r           = 0.0505;
% tau_d           = random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
tau_d           = (rand(numNeuron, 1)-0.5)*2*std_d*3 + median_d; %random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    paramsSet   = [Fm, Ca0, n, tau_d(nNeuron), tau_r, intNoise];
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model3(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 5) = ref_peak_list;

%% model 6 -- Nonlinearity
model_names{6} = 'Nonlinearity';
ref_spikeTimes_model3 = ref_spikeTimes;
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
randPairs       = randi([1 length(nlParaMat)], numNeuron, 1);
truncatedNormal = truncate(makedist('Normal'), -1.5, 1.5);
Fm              = 21.3240;
Ca0             = nlParaMat(randPairs, 1);
n               = nlParaMat(randPairs, 2);
tau_r           = 0.0505;
tau_d           = 1.7064;
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    paramsSet   = [Fm, Ca0(nNeuron), n(nNeuron), tau_d, tau_r, intNoise];
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model3(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 6) = ref_peak_list;

%% model 7 -- all parameters
model_names{7} = 'All params';
ref_spikeTimes_model3 = ref_spikeTimes;
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
randPairs       = randi([1 length(nlParaMat)], numNeuron, 1);
truncatedNormal = truncate(makedist('Normal'), -1.5, 1.5);
Fm              = 21.3240;
Ca0             = nlParaMat(randPairs, 1);
n               = nlParaMat(randPairs, 2);
std_d           = 0.4588;
median_d        = 1.7064;
tau_r           = 0.0505;
tau_d           = (rand(numNeuron, 1)-0.5)*2*std_d*3 + median_d; %random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    paramsSet   = [Fm, Ca0(nNeuron), n(nNeuron), tau_d(nNeuron), tau_r, intNoise];
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model3(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 7) = ref_peak_list;

%% model 8 -- all parameters + firing rate
model_names{8} = 'All params & firing rates';
ref_spikeTimes_model3 = ref_spikeTimes_model2;
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
randPairs       = randi([1 length(nlParaMat)], numNeuron, 1);
truncatedNormal = truncate(makedist('Normal'), -1.5, 1.5);
Fm              = 21.3240;
Ca0             = nlParaMat(randPairs, 1);
n               = nlParaMat(randPairs, 2);
std_d           = 0.4588;
median_d        = 1.7064;
tau_r           = 0.0505;
% tau_d           = random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
tau_d           = (rand(numNeuron, 1)-0.5)*2*std_d*3 + median_d; %random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
caData_model  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    paramsSet   = [Fm, Ca0(nNeuron), n(nNeuron), tau_d(nNeuron), tau_r, intNoise];
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model3(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model(nNeuron, :) = ave_/max(ave_);
end
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
ref_peak_list_model(:, 8) = ref_peak_list;

%% summary
% ref_peak_list_model_factor = mean(ref_peak_list_model(:,1));
% figure
% boxplot(ref_peak_list_model/ref_peak_list_model_factor)
% set(gca, 'XTickLabels', model_names)
% ylabel('Peakiness')
% box off
% set(gca, 'tickDir', 'out')
% setPrint(8, 6, 'hist_peakiness', 'pdf')


ref_peak_list_model_factor = mean(ref_peak_list_model(:,1));
figure
hold on
bar(1:8, mean(ref_peak_list_model)/ref_peak_list_model_factor, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
errorbar(1:8, mean(ref_peak_list_model)/ref_peak_list_model_factor, ...
    std(ref_peak_list_model/ref_peak_list_model_factor), 'k')
set(gca, 'TickDir', 'out')
xlim([0.5, 8.5])
ylim([0 1.05])
set(gca, 'XTickLabels', model_names, 'XTick', 1:8)
ylabel('Peakiness')
box off
set(gca, 'tickDir', 'out')
setPrint(8, 6, 'bar_peakiness', 'pdf')