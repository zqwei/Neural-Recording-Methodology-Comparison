addpath('../Func');
setDir;

numTrials        = 50;
load([TempDatDir 'Modeled_6s_AAV_S1.mat']);
s2cDataSet = nDataSet;
s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
% s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

load([TempDatDir 'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
caFiringRatesIpsi= caFiringRates;
caFiringRatesIpsi(:, 1, :, :) = caFiringRates(:, 2, :, :);
caFiringRatesIpsi(:, 2, :, :) = caFiringRates(:, 1, :, :);
caFiringRatesMean = squeeze(mean(caFiringRates, 4));
caFiringRatesIpsiMean = squeeze(mean(caFiringRatesIpsi, 4));
caFiringRatesMean    = [squeeze(caFiringRatesMean(:, 1, :)), squeeze(caFiringRatesMean(:, 2, :))];
caFiringRatesIpsiMean= [squeeze(caFiringRatesIpsiMean(:, 1, :)), squeeze(caFiringRatesIpsiMean(:, 2, :))];
numCaTime            = size(caFiringRatesMean, 2)/2;
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, 1:numCaTime)), squeeze(s2cFiringRatesMean(:, 2, 1:numCaTime))];

rhoMat                 = corr(caFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');
rhoMatIpsi             = corr(caFiringRatesIpsiMean', s2cFiringRatesMean', 'type', 'Spearman');
[maxRho, maxIndex]   = max(rhoMat, [], 1);
[maxRhoIpsi, maxIndexIpsi]   = max(rhoMatIpsi, [], 1);

nlParams = nan(size(s2cFiringRatesMean, 1), 4);

for nUnit = 1:size(s2cFiringRatesMean, 1)
    xdata = s2cFiringRates(nUnit, :, 1:numCaTime, :);
    xdata = sort(xdata, 4);
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRates(maxIndex(nUnit), :, :, :);
    else
        ydata = caFiringRatesIpsi(maxIndexIpsi(nUnit), :, :, :);
    end
    ydata = sort(ydata, 4);
    init_para = [min(ydata(:)), max(ydata(:)), median(xdata(:)), xdata(:)\ydata(:)];
    try
        param = lsqcurvefit(f,init_para,xdata(:),ydata(:),...
                [min(ydata(:)), max(ydata(:))/2, median(xdata(:))/3, 0], ...
                [0, max(ydata(:)), median(xdata(:))*3, 20]);
    catch
         param = nan(1, 4);
    end
        if param(4)<=0 || param(3)<0
            param = nan(1, 4);
        end
    nlParams(nUnit, :) = param;
end

save([TempDatDir 'FineTunedNLParamsS1.mat'], 'nlParams');