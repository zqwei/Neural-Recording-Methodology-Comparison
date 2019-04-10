%
% toy model
%
% =======================================
%
% Ziqiang Wei
% weiz@janelia.hhmi.org

% sample rate 
fs        = 1000;
t         = 1/fs:1/fs:4;
numNeuron = 100;
t_period  = 0.05; % rotation per t_period second
numSwitch = 4/t_period;
numTime   = t_period*fs;
numT      = numSwitch*numTime;

% coding direction is defined as act_no - act_yes
coding_direction = randn(numNeuron, numSwitch)*1;
figure,
imagesc(corr(coding_direction));
title('Toy model -- coding directions')

neural_act       = coding_direction - min(coding_direction, [], 2);
neural_act_yes   = ones(numNeuron, numSwitch)*0.5;
neural_act_no    = neural_act*0.4 + neural_act_yes;

% ephys data
numTrial = 100;
neural_act_yes = repmat(neural_act_yes, [1, 1, numTime]);
neural_act_yes = reshape(permute(neural_act_yes, [1, 3, 2]), numNeuron, []);
neural_act_no  = repmat(neural_act_no, [1, 1, numTime]);
neural_act_no  = reshape(permute(neural_act_no, [1, 3, 2]), numNeuron, []);

act_yes = zeros(numTrial, numNeuron, numT);
act_no  = act_yes;

for nTrial = 1:numTrial
    act_yes(nTrial, :, :) = neural_act_yes/fs > rand(size(neural_act_yes));
    act_no(nTrial, :, :)  = neural_act_no/fs > rand(size(neural_act_no));
end


last_period     = 0.2; % last period of decoding
last_period_bin = last_period * fs;

ave_yes         = squeeze(mean(act_yes(:, :, numT-last_period_bin:numT), 3));
ave_no          = squeeze(mean(act_no(:, :, numT-last_period_bin:numT), 3));
totTargets      = [ones(numTrial, 1); zeros(numTrial, 1)];
ave_            = [ave_yes; ave_no];

Mdl  = fitcdiscr(ave_, totTargets, 'discrimType','pseudoLinear', 'Kfold', 5);
performance = 1 - Mdl.kfoldLoss;
disp(['Ephys performance at last period of ' num2str(last_period) ' sec is ' num2str(performance)])


% S2C model
ca_linear   = zeros(numTrial*2, numNeuron, numT);
params = [1.7, 0.06, 3.0]; % decay, rise, internal noise


for nNeuron = 1:numNeuron
    spikeTimes = cell(numTrial*2, 1);
    for nTrial = 1:numTrial
        spikeTimes{nTrial}          = find(squeeze(act_yes(nTrial, nNeuron, :)))/fs;
        spikeTimes{nTrial+numTrial} = find(squeeze(act_no(nTrial, nNeuron, :)))/fs;
    end
    ca_linear(:, nNeuron, :) = spikeTimeToLinearImaging(spikeTimes, t, params);
end

g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
p_nonlinear = [0, 3, 4, 1];
ca_act = g(p_nonlinear, ca_linear) + randn(size(ca_linear)) * 0.3;


Mdl  = fitcdiscr(squeeze(ca_linear(:, :, end)), totTargets, 'discrimType','pseudoLinear', 'Kfold', 5);
performance = 1 - Mdl.kfoldLoss;
disp(['Ca++ performance (linear model) at last time is ' num2str(performance)])

Mdl  = fitcdiscr(squeeze(ca_act(:, :, end)), totTargets, 'discrimType','pseudoLinear', 'Kfold', 5);
performance = 1 - Mdl.kfoldLoss;
disp(['Ca++ performance (nonlinear model) at last time is ' num2str(performance)])
