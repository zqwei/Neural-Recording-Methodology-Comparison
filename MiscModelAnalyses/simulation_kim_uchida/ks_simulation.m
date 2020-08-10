addpath('../Func');
setDir;
TempDatDir = '../Backups/TempDat_2019_01_28/';

numNeuron  = 2000;
numTrial   = 1;
frameRate  = 14.84;
timeSeries = (-frameRate:frameRate*5)/frameRate;
intNoise     = 0.1;
extNoise     = 0.0;
fix_spk_num  = 2;
burst_frame  = floor(5*frameRate);
base_rate    = 3;
ratio        = 1/2;

ref_spikeTimes = cell(numNeuron, numTrial);
for nTrial   = 1:numTrial
    for nNeuron = 1:numNeuron
        spkTime = [];
        for n_ = 1:length(timeSeries)
            if n_ == burst_frame
                spk_num = fix_spk_num;
            else
                t_n_ = timeSeries(n_);
                rate_ = base_rate+ratio*t_n_*(t_n_>0)*(n_<burst_frame);
                spk_num = rate_/frameRate > rand(1);
            end
            if spk_num>0
                spkTime_ = sort(rand(spk_num, 1))/frameRate + timeSeries(n_);
                spkTime = [spkTime; spkTime_]; 
            end
        end
        ref_spikeTimes{nNeuron, nTrial} = spkTime;
    end
end

figure
subplot(3, 1, 1)
hold on
for nNeuron = 1:numNeuron
    spkTime = ref_spikeTimes{nNeuron, 1};
    plot(spkTime, ones(length(spkTime), 1)*nNeuron, '.k')
end
xlim([0, 5])
ylim([0, 200])
box off
xlabel('Time (sec)')
ylabel('Neuron Index')


subplot(3, 1, 2)
hold on
spkNum = zeros(1, length(timeSeries));
for nNeuron = 1:numNeuron
    spkTime = ref_spikeTimes{nNeuron, 1};
    spkNum_  = histcounts(spkTime, [timeSeries, timeSeries(end)+1/frameRate]);
    spkNum   = spkNum + spkNum_;
end
plot([timeSeries(2:end), timeSeries(end)+1/frameRate], spkNum/numNeuron*frameRate, '-k')
xlim([0, 5])
ylim([0, 30])
box off
xlabel('Time (sec)')
ylabel('Spikes (/sec)')



load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
Fm              = 21.3240;
Ca0             = median(nlParaMat(:, 1));
n               = median(nlParaMat(:, 2));
median_r        = 0.0505;
median_d        = 1.7064;
paramsSet       = [Fm, Ca0, n, median_d, median_r, intNoise];
numTimeBin      = length(timeSeries);
caData          = nan(numNeuron, numTimeBin);
for nNeuron     = 1:numNeuron
    caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes(nNeuron, :), timeSeries, paramsSet, 3);
%     ave_        = mean(caNeuron, 1)+extNoise*randn(1, numTimeBin);
    f0          = mean(caNeuron(1:25));
    caData(nNeuron, :) = caNeuron/f0 -1; %/max(caNeuron)
end

subplot(3, 1, 3)
plot(timeSeries, mean(caData, 1), '-k')
xlim([0, 5])
ylim([0, 0.75])
box off
xlabel('Time (sec)')
ylabel('average dF/F')

setPrint(8, 12, ['simulation_kim_uchida_v' num2str(fix_spk_num)], 'pdf')
