addpath('../Func');
setDir;

numTrials        = 50;
load([TempDatDir 'Modeled_6s_AAV_Intra.mat']);
s2cDataSet = nDataSet;
s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

load([TempDatDir 'Shuffle_Ca_Slow_Short_Delay_Virus_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
caFiringRatesIpsi= caFiringRates;
caFiringRatesIpsi(:, 1, :, :) = caFiringRates(:, 2, :, :);
caFiringRatesIpsi(:, 2, :, :) = caFiringRates(:, 1, :, :);
caFiringRatesMean = squeeze(mean(caFiringRates, 4));
caFiringRatesIpsiMean = squeeze(mean(caFiringRatesIpsi, 4));
caFiringRatesMean    = [squeeze(caFiringRatesMean(:, 1, :)), squeeze(caFiringRatesMean(:, 2, :))];
caFiringRatesIpsiMean= [squeeze(caFiringRatesIpsiMean(:, 1, :)), squeeze(caFiringRatesIpsiMean(:, 2, :))];

rhoMat                 = corr(caFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');
rhoMatIpsi             = corr(caFiringRatesIpsiMean', s2cFiringRatesMean', 'type', 'Spearman');
[maxRho, maxIndex]   = max(rhoMat, [], 1);
[maxRhoIpsi, maxIndexIpsi]   = max(rhoMatIpsi, [], 1);

f = @(p,x) p(1) + p(2)./ (1 + exp((p(3)-x)*p(4)));
nlParams = nan(2, size(s2cFiringRatesMean, 1), 4);
for nUnit = 1:size(s2cFiringRatesMean, 1)
    xdata = s2cFiringRates(nUnit, :, :, :);
    xdata = sort(xdata, 4);
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRates(maxIndex(nUnit), :, :, :);
    else
        ydata = caFiringRatesIpsi(maxIndexIpsi(nUnit), :, :, :);
    end
    ydata = sort(ydata, 4);
    init_para = [min(ydata(:)), max(ydata(:)), median(xdata(:)), xdata(:)\ydata(:)];
    % param = nlinfit(xdata(:),ydata(:),f,init_para);
    param = lsqcurvefit(f,init_para,xdata(:),ydata(:),...
            [min(ydata(:)), max(ydata(:))/2, median(xdata(:))/3, 0], ...
            [0, max(ydata(:)), median(xdata(:))*3, 20]);
    % [param,~]=sigm_fit(xdata(:), ydata(:), [], [], false);
    % param(2) = param(2) - param(1);
    % param(4) = param(4) * log(10);
    if param(4)<=0 || param(3)<0
        param = nan(1, 4);
    end
    nlParams(1, nUnit, :) = param;
%     nls2cFiringRates(nUnit, :, :, :) = f(param, s2cFiringRates(nUnit, :, :, :));
end

numTrials        = 50;
load([TempDatDir 'Modeled_GP43_Intra.mat']);
s2cDataSet = nDataSet;
s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

load([TempDatDir 'Shuffle_Ca_Slow_Short_Delay_withOLRemoval.mat']);
caDataSet        = nDataSet;
caFiringRates    = generateDPCAData(caDataSet, numTrials);
caFiringRatesIpsi= caFiringRates;
caFiringRatesIpsi(:, 1, :, :) = caFiringRates(:, 2, :, :);
caFiringRatesIpsi(:, 2, :, :) = caFiringRates(:, 1, :, :);
caFiringRatesMean = squeeze(mean(caFiringRates, 4));
caFiringRatesIpsiMean = squeeze(mean(caFiringRatesIpsi, 4));
caFiringRatesMean    = [squeeze(caFiringRatesMean(:, 1, :)), squeeze(caFiringRatesMean(:, 2, :))];
caFiringRatesIpsiMean= [squeeze(caFiringRatesIpsiMean(:, 1, :)), squeeze(caFiringRatesIpsiMean(:, 2, :))];

rhoMat                 = corr(caFiringRatesMean', s2cFiringRatesMean', 'type', 'Spearman');
rhoMatIpsi             = corr(caFiringRatesIpsiMean', s2cFiringRatesMean', 'type', 'Spearman');
[maxRho, maxIndex]   = max(rhoMat, [], 1);
[maxRhoIpsi, maxIndexIpsi]   = max(rhoMatIpsi, [], 1);

for nUnit = 1:size(s2cFiringRatesMean, 1)
    xdata = s2cFiringRates(nUnit, :, :, :);
    xdata = sort(xdata, 4);
    if maxRho(nUnit)>maxRhoIpsi(nUnit)
        ydata = caFiringRates(maxIndex(nUnit), :, :, :);
    else
        ydata = caFiringRatesIpsi(maxIndexIpsi(nUnit), :, :, :);
    end
    ydata = sort(ydata, 4);
    init_para = [min(ydata(:)), max(ydata(:)), median(xdata(:)), xdata(:)\ydata(:)];
    param = lsqcurvefit(f,init_para,xdata(:),ydata(:),...
            [min(ydata(:)), max(ydata(:))/2, median(xdata(:))/3, 0], ...
            [0, max(ydata(:)), median(xdata(:))*3, 20]);
    if param(4)<=0 || param(3)<0
        param = nan(1, 4);
    end
    nlParams(2, nUnit, :) = param;
end

save([TempDatDir 'FineTunedNLParamsIntra.mat'], 'nlParams');