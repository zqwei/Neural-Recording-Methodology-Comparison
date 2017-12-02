
function exampleRampingDownNeuron
    addpath('../Func');
    setDir;

    if ~exist([PlotDir 'SingleUnitsRampingDown'],'dir')
        mkdir([PlotDir 'SingleUnitsRampingDown'])
    end    
    
    
    load ([TempDatDir 'DataListShuffle.mat']);
    nData = 1; % plot raster and psth
    load([TempDatDir DataSetList(nData).name '.mat'])
    spikeDataSet = nDataSet;   
    params = DataSetList(nData).params;
    
    yesRampDown          = zeros(length(nDataSet), 1);
    noRampDown           = zeros(length(nDataSet), 1);
    
    timePoints(1)        = sum(DataSetList(nData).params.timeSeries<DataSetList(nData).params.polein);
    timePoints(2)        = sum(DataSetList(nData).params.timeSeries<0);
    
    for nUnit            = 1:length(spikeDataSet)
        meanPreSample    = [mean(spikeDataSet(nUnit).unit_yes_trial(:, 1:timePoints(1)), 2); mean(spikeDataSet(nUnit).unit_no_trial(:, 1:timePoints(1)), 2)];
        meanYesSample    = mean(spikeDataSet(nUnit).unit_yes_trial(:, timePoints(1):timePoints(2)), 2);
        meanNoSample     = mean(spikeDataSet(nUnit).unit_no_trial(:, timePoints(1):timePoints(2)), 2);
        if ttest2(meanPreSample, meanYesSample, 'tail', 'right')
            yesRampDown(nUnit) = 1;
        elseif ttest2(meanPreSample, meanYesSample, 'tail', 'left')
            yesRampDown(nUnit) = -1;
        end
        
        if ttest2(meanPreSample, meanNoSample, 'tail', 'right')
            noRampDown(nUnit)  = 1;
        elseif ttest2(meanPreSample, meanNoSample, 'tail', 'left')
            noRampDown(nUnit)  = -1;
        end
    end
    

    nData = 4; % plot ca comparison
    load([TempDatDir DataSetList(nData).name '.mat'])
    caDataSet    = nDataSet;
    load ([TempDatDir 'DataListS2CModel.mat']);    
    nData = 4; % plot s2c model
    load([TempDatDir DataSetList(nData).name '.mat'])
    s2cDataSet   = nDataSet;
    


    dynamicalNeuronIndex = find(yesRampDown==1 & noRampDown==1);
    mCellSet             = find(yesRampDown==1 & noRampDown==1);
    
    for nCellid = 1:length(dynamicalNeuronIndex)
        nCell   = dynamicalNeuronIndex(nCellid);
        figure;
        % spike
        subplot(2, 3, 1)
        plotRaster(spikeDataSet(nCell), params);
        subplot(2, 3, 4)
        plotPSTH(spikeDataSet(nCell), params, 'Firing rate (Hz)');
        
        %
        subplot(2, 3, 3)
        plotDff(s2cDataSet(nCell), params)
        subplot(2, 3, 6)
        plotPSTH(s2cDataSet(nCell), params, 'DF/F');
        
%         mCell = findSimilarCellToS2CModel(meanS2CDataSet(nCell,:), meanCaDataSet);
        mCell = mCellSet(nCellid);
        subplot(2, 3, 2)
        plotDff(caDataSet(mCell), params)
        subplot(2, 3, 5)
        plotPSTH(caDataSet(mCell), params, 'DF/F');
        setPrint(8*3, 6*2, [PlotDir 'SingleUnitsRampingDown/SingleUnitsRampingDownExampleNeuron_' num2str(nCell, '%04d')])
        close all
    end

end

function plotRaster(spikeDataSet, params)
    spkTimes{1}    = spikeDataSet.unit_yes_trial_spk_time;
    spkTimes{2}    = spikeDataSet.unit_no_trial_spk_time;
    hold on;
    color_index    = [0.7  0 0; 0 0 0.7];
    spkLoc         = 0;
    for nPlot            = 1:2
        hold on;
        for nTrial       = 1:20
            spkLoc       = spkLoc + 1;
            spikeTimes   = spkTimes{nPlot}{nTrial};
            spikeTrial   = ones(length(spikeTimes), 1) * spkLoc;
            plot(spikeTimes, spikeTrial, '.', 'color', color_index(nPlot, :));
        end
        spkLoc     = spkLoc + 1;%length(spkTimes{nPlot}) + 3;
    end

    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 41])
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end

function plotDff(spikeDataSet, params)
    sigma                         = 0.15 / params.binsize; % 200 ms
    filterLength                  = 11;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse); 

    nUnitData        = spikeDataSet.unit_yes_trial;
    yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
    nUnitData        = spikeDataSet.unit_no_trial;
    noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
    cmin             = min([mean(yesUnitData), mean(noUnitData)]);
    cmax             = max([mean(yesUnitData), mean(noUnitData)]);
    


    hold on
    actMat = nan(41, size(spikeDataSet.unit_yes_trial, 2));
    actMat(1:20, :) = yesUnitData(1:20, :);
    actMat(22:end, :) = noUnitData(1:20, :);
    imagesc(params.timeSeries, 1:41, actMat);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','w','Linestyle','--','linewid', 1.0);
    hold off;
    ylim([1 41])
    colormap(gray)
    
    caxis([cmin, cmax]);
%     colorbar
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    axis off
end


function plotPSTH(spikeDataSet, params, ylabelName)
    sigma                         = 0.15 / params.binsize; % 200 ms
    filterLength                  = 11;
    filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
    filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
    filterInUse                   = filterInUse / sum (filterInUse); 

    color_index    = [0.7  0 0; 0 0 0.7];
    hold on;
    nUnitData        = spikeDataSet.unit_yes_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(1, :)}, 0.5);
    nUnitData        = spikeDataSet.unit_no_trial;
    nUnitData        = getGaussianPSTH (filterInUse, nUnitData, 2);
    shadedErrorBar(params.timeSeries, mean(nUnitData, 1),...
        std(nUnitData, 1)/sqrt(size(nUnitData, 1)),...
        {'k-', 'linewid', 1.0, 'color', color_index(2, :)}, 0.5);
    gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
    hold off;
    ylabel(ylabelName);
    xlabel('Time (s)');
    xlim([params.timeSeries(1) params.timeSeries(end)]);
    set(gca, 'TickDir', 'out')
end


function meanDataSet = meanData(nDataSet)
    meanDataSet = nan(length(nDataSet), size(nDataSet(1).unit_yes_trial, 2)*2);
    for nCell = 1:length(nDataSet)
        meanDataSet(nCell, :) = [mean(nDataSet(nCell).unit_yes_trial), mean(nDataSet(nCell).unit_no_trial)];
    end
end

function mIndex = findSimilarCellToS2CModel(nCellS2CDataSet, caDataSet)
    corrDat     = corr(nCellS2CDataSet', caDataSet');
    [~, mIndex] = max(corrDat);
    mIndex      = mIndex(1);
end