addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

nData     = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
neuronRemoveList = false(length(nDataSet), 1);

disp(DataSetList(nData).name)
params      = DataSetList(nData).params;
timeBins    = params.timeSeries;
binsize     = params.binsize;
timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
noProfileMatrix     = yesProfileMatrix;
positivePeak        = false(length(nDataSet),1);

for nUnit        = 1:length(nDataSet)
    yesData      = mean(nDataSet(nUnit).unit_yes_trial);
    noData       = mean(nDataSet(nUnit).unit_no_trial);
    maxData      = max([yesData, noData]);
    minData      = min([yesData, noData]);
    rData        = (maxData - minData);
    yesData      = (yesData - minData)/(maxData - minData);
    noData       = (noData - minData)/(maxData - minData);
    yesProfileMatrix(nUnit, :)    = yesData;
    noProfileMatrix(nUnit, :)     = noData;
    positivePeak(nUnit)       = mean(yesData(timePoints(1):timePoints(2))) <= mean(yesData(timePoints(2):timePoints(4))) ...
                               || mean(noData(timePoints(1):timePoints(2))) <= mean(noData(timePoints(2):timePoints(4)));

end
actMat        = [yesProfileMatrix, noProfileMatrix];
actMat        = actMat(positivePeak, :);
[~, maxId]    = max(actMat, [], 2);


peakOnlyData  = nDataSet(positivePeak);
for nCell     = 1:length(peakOnlyData)
    time_     = maxId(nCell);
    numYes    = size(peakOnlyData(nCell).unit_yes_trial, 1);
    numNo     = size(peakOnlyData(nCell).unit_no_trial, 1);
    if time_  <= numTimeBin
        peak_time = timeBins(time_);
        spk_time  = peakOnlyData(nCell).unit_yes_trial_spk_time;
        for ntrial = 1:numYes
            spk_trial = spk_time{ntrial};
            spk_trial = spk_trial(spk_trial>peak_time-binsize & spk_trial<peak_time+binsize);
            spk_time{ntrial} = spk_trial;
            spk_count = zeros(1, numTimeBin);
            spk_count(time_) = length(spk_trial);
            peakOnlyData(nCell).unit_yes_trial(ntrial, :) = spk_count;
        end
        peakOnlyData(nCell).unit_yes_trial_spk_time = spk_time;
        peakOnlyData(nCell).unit_no_trial = zeros(numNo, numTimeBin);
        peakOnlyData(nCell).unit_no_trial_spk_time = cell(numNo, 1);
    else
        time_     = time_ - numTimeBin;
        peak_time = timeBins(time_);
        spk_time  = peakOnlyData(nCell).unit_no_trial_spk_time;
        for ntrial = 1:numNo
            spk_trial = spk_time{ntrial};
            spk_trial = spk_trial(spk_trial>peak_time-binsize & spk_trial<peak_time+binsize);
            spk_time{ntrial} = spk_trial;
            spk_count = zeros(1, numTimeBin);
            spk_count(time_) = length(spk_trial);
            peakOnlyData(nCell).unit_no_trial(ntrial, :) = spk_count;
        end
        peakOnlyData(nCell).unit_no_trial_spk_time = spk_time;
        peakOnlyData(nCell).unit_yes_trial = zeros(numYes, numTimeBin);
        peakOnlyData(nCell).unit_yes_trial_spk_time = cell(numYes, 1);
    end
end

yesProfileMatrix    = nan(length(peakOnlyData), numTimeBin);
noProfileMatrix     = yesProfileMatrix;

for nUnit        = 1:length(peakOnlyData)
    yesData      = mean(peakOnlyData(nUnit).unit_yes_trial);
    noData       = mean(peakOnlyData(nUnit).unit_no_trial);
    maxData      = max([yesData, noData]);
    minData      = min([yesData, noData]);
    rData        = (maxData - minData);
    yesData      = (yesData - minData)/(maxData - minData);
    noData       = (noData - minData)/(maxData - minData);
    yesProfileMatrix(nUnit, :)    = yesData;
    noProfileMatrix(nUnit, :)     = noData;
end
actMat        = [yesProfileMatrix, noProfileMatrix];
[~, maxId]    = max(actMat, [], 2);
[peakId, sortId]   = sort(maxId);
% figure;
% imagesc(actMat(sortId, :))
figure;
plot(peakId, 1:length(peakOnlyData), 'o')

load([TempDatDir 'nlParaMat.mat'],'nlParaMat');
load([TempDatDir 'FineTunedNLParams.mat'], 'nlParams');
truncatedNormal        = truncate(makedist('Normal'), -1.5, 1.5);
params.Fm              = random(truncatedNormal, length(nDataSet), 1) * 10.6673 + 21.3240;
randPairs              = randi([1 length(nlParaMat)], length(nDataSet), 1);
params.Ca0             = nlParaMat(randPairs, 1);
params.n               = nlParaMat(randPairs, 2);
std_r                  = 0.0246;
median_r               = 0.0505;
std_d                  = 0.4588;
median_d               = 1.7064;
params.tau_r           = random(truncatedNormal, length(nDataSet), 1) *  std_r + median_r;
params.tau_d           = random(truncatedNormal, length(nDataSet), 1) *  std_d + median_d;
params.intNoise        = 1.5;
params.extNoise        = 1.5;
nDataSet               = getFakeCaImagingDataSigmoid(peakOnlyData, params);


% nonlinear function
g = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));

int_noise = 2.0;
ext_noise = 0.0;

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

for nUnit        = 1:length(nDataSet)
    yesData      = mean(nDataSet(nUnit).unit_yes_trial);
    noData       = mean(nDataSet(nUnit).unit_no_trial);
    maxData      = max([yesData, noData]);
    minData      = min([yesData, noData]);
    rData        = (maxData - minData);
    yesData      = (yesData - minData)/(maxData - minData);
    noData       = (noData - minData)/(maxData - minData);
    yesProfileMatrix(nUnit, :)    = yesData;
    noProfileMatrix(nUnit, :)     = noData;
end
actMat        = [yesProfileMatrix, noProfileMatrix];
[~, maxId]    = max(actMat, [], 2);
[peakId, sortId]   = sort(maxId);
% figure;
% imagesc(actMat(sortId, :))
% figure;
hold on
plot(peakId, 1:length(nDataSet), 'o')
hold off
box off
set(gca, 'tickdir', 'out');
xlabel('Time')
ylabel('Neuron index');
box off
legend('Ephys', 'S2C')
setPrint(8, 6, 'peakLocationSinglePeakS2C', 'pdf')


timeStep  = DataSetList(nData).params.timeSeries;
timeTag   = timePoints(2):timePoints(4)+13; % sample to response
numTime   = length(timeTag);
polein    = DataSetList(nData).params.polein;
poleout   = DataSetList(nData).params.poleout;

countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
std(countMaxId([timeTag, timeTag+numTimeBin]))%/mean(countMaxId([timeTag, timeTag+77]))
[bootstat,bootsam] = bootstrp(1000,@std,countMaxId([timeTag, timeTag+numTimeBin]));
std(bootstat)
figure;
hold on;
bplot = bar(timeStep, countMaxId(1:numTimeBin), 1, 'facecolor', 'b', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
bplot = bar(timeStep, countMaxId(1+numTimeBin:end), 1, 'facecolor', 'r', 'edgecolor', 'none');
bplot.FaceAlpha = 0.5;
xlim([timeStep(1) params.timeSeries(end)])
ylim([0 8])
gridxy ([polein, poleout, 0],[1/(numTimeBin-8)/2*100], 'Color','k','Linestyle','--','linewid', 1.0)
box off
set(gca,'TickDir','out')
ylabel('% Max peak')
xlabel('Time')
hold off
setPrint(8, 6, 'peakLocationSinglePeakS2C_distrs', 'pdf')