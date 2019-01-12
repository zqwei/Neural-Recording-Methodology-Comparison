addpath('../Func');
setDir;

numNeuron  = 50;
numTrial   = 30;
frameRate  = 14.84;
numTimeBin = 61;
timeSeries = (0:60)/frameRate;
peakTime   = zeros(1, numTimeBin);
peakTime(11) = 1;
ephys      = sqrt(mean((peakTime - 1/numTimeBin).^2))*100;

ref_spikeTimes = cell(numNeuron, numTrial);
for nTrial   = 1:numTrial
    spkTime  = sort(rand(3, 1))/frameRate + timeSeries(10);
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
intNoise        = 0.3;
extNoise        = 0.0;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData          = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData(nNeuron, :) = ave_/max(ave_);
end

[~, maxId] = max(caData, [], 2);
countMaxId = hist(maxId, 1:numTimeBin);
ref_peak = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
% std based on jackknife
ref_peak_list = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
figure;
imagesc(caData)
title('Control')

% model one, spike time varies
ref_spikeTimes_model1 = cell(numNeuron, numTrial);
for nTrial   = 1:numTrial
    for nNeuron = 1:numNeuron
        spkTime  = sort(rand(3, 1))/frameRate + timeSeries(10);
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
intNoise        = 0.3;
extNoise        = 0.0;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData_model1  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model1(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model1(nNeuron, :) = ave_/max(ave_);
end

[~, maxId] = max(caData_model1, [], 2);
countMaxId = hist(maxId, 1:numTimeBin);
ref_peak_model1 = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
% std based on jackknife
ref_peak_list_model1 = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model1;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list_model1(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
figure;
imagesc(caData)
title('Spike times vs Control')

% model two, firing rate varies
ref_spikeTimes_model2 = cell(numNeuron, numTrial);
rate_        = 0.1:0.1:20;
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
intNoise        = 0.3;
extNoise        = 0.0;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
caData_model2  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model2(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model2(nNeuron, :) = ave_/max(ave_);
end

[~, maxId] = max(caData_model2, [], 2);
countMaxId = hist(maxId, 1:numTimeBin);
ref_peak_model2 = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
% std based on jackknife
ref_peak_list_model2 = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model2;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list_model2(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
figure;
imagesc(caData_model2)
title('Firing rate vs Control')


% model three, firing rate varies
ref_spikeTimes_model3 = cell(numNeuron, numTrial);
load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
randPairs       = randi([1 length(nlParaMat)], numNeuron, 1);
truncatedNormal = truncate(makedist('Normal'), -1.5, 1.5);
Fm              = 21.3240;
Ca0             = nlParaMat(randPairs, 1);
n               = nlParaMat(randPairs, 2);
std_r           = 0.0246;
median_r        = 0.0505;
std_d           = 0.4588;
median_d        = 1.7064;
tau_r           = random(truncatedNormal, numNeuron, 1) *  std_r + median_r;
tau_d           = random(truncatedNormal, numNeuron, 1) *  std_d + median_d;
intNoise        = 0.3;
extNoise        = 0.0;
caData_model3  = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    paramsSet   = [Fm, Ca0(nNeuron), n(nNeuron), tau_d(nNeuron), tau_r(nNeuron), intNoise];
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes_model2(nNeuron, :), timeSeries, paramsSet, 0);
    ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    ave_(ave_<0)=0;
    caData_model3(nNeuron, :) = ave_/max(ave_);
end

[~, maxId] = max(caData_model3, [], 2);
countMaxId = hist(maxId, 1:numTimeBin);
ref_peak_model3 = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
% std based on jackknife
ref_peak_list_model3 = nan(numNeuron, 1);
for nNeuron = 1:numNeuron
    ca_    = caData_model2;
    ca_(nNeuron, :) = [];
    [~, maxId] = max(ca_, [], 2);
    countMaxId = hist(maxId, 1:numTimeBin);
    ref_peak_list_model3(nNeuron) = sqrt(mean((countMaxId - 1).^2))/numTimeBin*100;
end
figure;
imagesc(caData_model3)
title('S2C params vs Control')

figure
boxplot([ref_peak_list, ref_peak_list_model1, ref_peak_list_model2, ref_peak_list_model3])