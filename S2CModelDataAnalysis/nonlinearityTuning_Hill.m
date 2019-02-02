addpath('../Func');
setDir;

%%% define nonlinearity
f = @(p,x) p(1) + p(2)*x.^p(3)./(p(4)+x.^p(3));

%%% Datasets to fit
sourceDataSetName = {'Modeled_6s_AAV', ...
                     'Modeled_GP43', ...
                     'Modeled_GP517', ...
                     'Modeled_ALM_Long_GP43', ...
                     'Modeled_ALM_Long_GP517', ...
                     'Modeled_6s_AAV_S1'};
targetDataSetName = {'Shuffle_Ca_Slow_Short_Delay_Virus', ...
                     'Shuffle_Ca_Slow_Short_Delay', ...
                     'Shuffle_Ca_Fast_SShort_Delay', ...
                     'Shuffle_Ca_Slow_Long_Delay', ...
                     'Shuffle_Ca_Fast_Long_Delay', ...
                     'Shuffle_Ca_SP_S1_Slow_Short_Delay_Virus'};

maxNumNeuron = 4000;
numParam = 4;
nlParams = nan(length(sourceDataSetName), maxNumNeuron, numParam);           

for nData        = 1:length(sourceDataSetName)
    numTrials        = 50;
    load([TempDatDir sourceDataSetName{nData} '.mat']);
    s2cDataSet = nDataSet;
    s2cFiringRates   = generateDPCADataLinearData(s2cDataSet, numTrials);
    s2cFiringRatesMean   = squeeze(mean(s2cFiringRates, 4));
%     s2cFiringRatesMean   = [squeeze(s2cFiringRatesMean(:, 1, :)), squeeze(s2cFiringRatesMean(:, 2, :))];

    load([TempDatDir targetDataSetName{nData} '_withOLRemoval.mat']);
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

    for nUnit = 1:size(s2cFiringRatesMean, 1)
        xdata = s2cFiringRates(nUnit, :, 1:numCaTime, :);
        xdata = sort(xdata, 4);
        if maxRho(nUnit)>maxRhoIpsi(nUnit)
            ydata = caFiringRates(maxIndex(nUnit), :, :, :);
        else
            ydata = caFiringRatesIpsi(maxIndexIpsi(nUnit), :, :, :);
        end
        ydata = sort(ydata, 4);
        init_para = [min(ydata(:)), max(ydata(:)), 2, median(ydata(:))];
        try
            param = lsqcurvefit(f,init_para,xdata(:),ydata(:),...
                    [min(ydata(:)), max(ydata(:))/2, 1, 0], ...
                    [0, max(ydata(:)), 4, inf]);
        catch
             param = nan(1, numParam);
        end
        if param(4)<=0 || param(3)<0
            param = nan(1, 4);
        end
        nlParams(nData, nUnit, :) = param;
    end
end


save([TempDatDir 'FineTunedNLParams_Hill.mat'], 'nlParams');