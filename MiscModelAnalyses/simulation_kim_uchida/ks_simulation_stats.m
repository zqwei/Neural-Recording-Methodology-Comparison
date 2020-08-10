addpath('../Func');
setDir;
TempDatDir = '../Backups/TempDat_2019_01_28/';

numNeuron  = 2000;
numTrial   = 1;
frameRate  = 14.84;
timeSeries = (-frameRate:frameRate*5)/frameRate;
intNoise     = 0.1;
extNoise     = 0.0;
burst_frame  = floor(5*frameRate);

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


base_rate_list   = [1 2 3 4 8 10];
ratio_list       = 0.5:0.5:3;
fix_spk_num_list = 1:5;

% spk_rate_list = zeros(length(base_rate_list), length(ratio_list), length(fix_spk_num_list), numTimeBin);
% ca_list = zeros(length(base_rate_list), length(ratio_list), length(fix_spk_num_list), numTimeBin);
% 
% for nbase = 1:length(base_rate_list)
%     for nratio = 1:length(ratio_list)
%         for nfix = 1:length(fix_spk_num_list)
%             base_rate   = base_rate_list(nbase);
%             ratio  = ratio_list(nratio);
%             fix_spk_num = fix_spk_num_list(nfix);
%             ref_spikeTimes = cell(numNeuron, numTrial);
%             for nTrial   = 1:numTrial
%                 for nNeuron = 1:numNeuron
%                     spkTime = [];
%                     for n_ = 1:length(timeSeries)
%                         if n_ == burst_frame
%                             spk_num = fix_spk_num;
%                         else
%                             t_n_ = timeSeries(n_);
%                             rate_ = base_rate+ratio*t_n_*(t_n_>0)*(n_<burst_frame);
%                             spk_num = rate_/frameRate > rand(1);
%                         end
%                         if spk_num>0
%                             spkTime_ = sort(rand(spk_num, 1))/frameRate + timeSeries(n_);
%                             spkTime = [spkTime; spkTime_]; 
%                         end
%                     end
%                     ref_spikeTimes{nNeuron, nTrial} = spkTime;
%                 end
%             end
% 
%             spkNum = zeros(1, length(timeSeries));
%             for nNeuron = 1:numNeuron
%                 spkTime = ref_spikeTimes{nNeuron, 1};
%                 spkNum_  = histcounts(spkTime, [timeSeries, timeSeries(end)+1/frameRate]);
%                 spkNum   = spkNum + spkNum_;
%             end
%             spk_rate = spkNum/numNeuron*frameRate;
%             spk_rate_list(nbase, nratio, nfix, :)=spk_rate;
% 
%             for nNeuron     = 1:numNeuron
%                 caNeuron    = spikeTimeToImagingSigmoid(ref_spikeTimes(nNeuron, :), timeSeries, paramsSet, 3);
%                 f0          = mean(caNeuron(1:25));
%                 caData(nNeuron, :) = caNeuron/f0 -1; %/max(caNeuron)
%             end
%             ca_ = mean(caData, 1);
%             ca_list(nbase, nratio, nfix, :)=ca_;
%         end
%     end
% end
% 
% save('tempDat', 'spk_rate_list', 'ca_list')


load('tempDat')
spk_rp_ratio = max(spk_rate_list, [], 4)./max(spk_rate_list(:, :, :, 1:burst_frame-1), [], 4);
ca_rp_ratio = max(ca_list, [], 4)./max(ca_list(:, :, :, 1:burst_frame-1), [], 4);

figure;
hold on
for n_base = 1:length(base_rate_list)
    scatter(squeeze(spk_rp_ratio(n_base,1,:))-1, squeeze(ca_rp_ratio(n_base,1,:))-1, fix_spk_num_list*20+6, 'filled')
end
plot([1e-5, 100], [1e-5, 100], '--k')
xlim([0.1, 100])
set(gca,'xscale','log', 'yscale', 'log')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4.5 3]);
saveas(gcf, 'spk_ca_base_fix_spk.svg')

figure;
hold on
for n_base = 1:length(base_rate_list)
    scatter(mean(spk_rp_ratio(n_base, :, :), 3)-1, mean(ca_rp_ratio(n_base, :, :), 3)-1, ratio_list*10+3, 'filled')
end
% plot([1e-5, 100], [1e-5, 100], '--k')
xlim([0.5, 20])
set(gca,'xscale','log', 'yscale', 'log')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4.5 3]);
saveas(gcf, 'spk_ca_base_ramp_speed.svg')

