%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPCA and PCA across trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

combinedParams = {{1}, {2}, {[1 2]}};
margNames      = {'Stim', 'Time', 'Inter'};
if ~exist([PlotDir 'CollectedUnitsPCA'],'dir')
    mkdir([PlotDir 'CollectedUnitsPCA'])
end


numComps       = 10;
cmap = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

ROCThres       = 0.50;
nDataList      = [1 3 4 13 10 11 12];
evMat          = zeros(length(nDataList), length(combinedParams), numComps);
numTrials      = 100;

for mData      = 1:length(nDataList)
    nData      = nDataList(mData);
    if nData == 13
        load([TempDatDir DataSetList(1).name '.mat'])
        selectedNeuronalIndex = DataSetList(1).ActiveNeuronIndex';
        depth_list = [nDataSet.depth_in_um];
        selectedNeuronalIndex = selectedNeuronalIndex & depth_list < 471;
        nData = 1;
    else
        load([TempDatDir DataSetList(nData).name '.mat'])
        selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    end
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    totalVar           = sum(sum(pcaX.^2));
    [~, S, Wpca]       = svd(pcaX');

    PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
    for i=1:length(Xmargs)
        PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
    end
    evMat(mData, :, :) = PCAmargVar(:, 1:numComps);
end

figure
plot(1:numComps, squeeze(cumsum(sum(evMat, 2), 3))', '-o')
legend({'Ephys', '6s-TG', '6s-AAV', 'Ephys', '6f-TG', 'S1', '6s-AAV'})
box off
ylim([0 1])
xlim([0.5 numComps+0.5]);
xlabel('PC components')
ylabel('EV')  
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'PCA_comps', 'pdf') 

figure
plot(1:numComps, squeeze(cumsum(sum(evMat, 2), 3))', '-o')
box off
ylim([0 1])
xlim([0.5 numComps+0.5]);
xlabel('PC components')
ylabel('EV')  
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'PCA_comps_no_legend', 'pdf') 


figure
for n_ = 1:3
    bar(1:length(nDataList), squeeze(evMat(:, :, n_)),'stacked', 'edgecolor', 'none');
    ax = gca;
    ax.XTick = 1:length(nDataList);
    ax.XTickLabel = {'Ephys', '6s-TG', '6s-AAV', 'Ephys', '6f-TG', 'S1', '6s-AAV'};
    ax.XTickLabelRotation = 30;
    ax.YTick = 0:0.1:0.7;
    box off
    ylim([0 0.7])
    xlim([0.5 length(nDataList)+0.5]);
    xlabel('PC components')
    ylabel('EV')  
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, ['PCA_n_comp_' num2str(n_)], 'pdf') 
end



load ([TempDatDir 'DataListS2CShuffle.mat']);
ROCThres       = 0.50;
nDataList      = [3 6 15 18];
evMat          = zeros(length(nDataList), length(combinedParams), numComps);
numTrials      = 100;

for mData      = 1:length(nDataList)
    nData      = nDataList(mData);
    load([TempDatDir DataSetList(nData).name '.mat'])
    selectedNeuronalIndex = DataSetList(nData).ActiveNeuronIndex';
    
    oldDataSet          = nDataSet;
    selectedNeuronalIndex = selectedHighROCneurons(oldDataSet, DataSetList(nData).params, ROCThres, selectedNeuronalIndex);
    nDataSet              = oldDataSet(selectedNeuronalIndex);
    
    firingRates        = generateDPCAData(nDataSet, numTrials);
    firingRatesAverage = nanmean(firingRates, ndims(firingRates));
    pcaX               = firingRatesAverage(:,:);
    firingRatesAverage = bsxfun(@minus, firingRatesAverage, mean(pcaX,2));
    pcaX               = bsxfun(@minus, pcaX, mean(pcaX,2));
    Xmargs             = dpca_marginalize(firingRatesAverage, 'combinedParams', combinedParams, 'ifFlat', 'yes');
    totalVar           = sum(sum(pcaX.^2));
    [~, S, Wpca]       = svd(pcaX');

    PCAmargVar         = zeros(length(combinedParams), length(nDataSet));
    for i=1:length(Xmargs)
        PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / totalVar;
    end
    evMat(mData, :, :) = PCAmargVar(:, 1:numComps);
end

figure
plot(1:numComps, squeeze(cumsum(sum(evMat, 2), 3))', '-o')
legend({'S2C 6s-TG', 'S2C 6s-AAV', 'S2C 6f-TG', 'S2C S1 6s-AAV'})
box off
ylim([0 1])
xlim([0.5 numComps+0.5]);
xlabel('PC components')
ylabel('EV')  
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'S2C_PCA_comps', 'pdf') 

figure
plot(1:numComps, squeeze(cumsum(sum(evMat, 2), 3))', '-o')
box off
ylim([0 1])
xlim([0.5 numComps+0.5]);
xlabel('PC components')
ylabel('EV')  
set(gca, 'TickDir', 'out')
setPrint(8, 6, 'S2C_PCA_comps_no_legend', 'pdf') 

figure
for n_ = 1:3
    bar(1:length(nDataList), squeeze(evMat(:, :, n_)),'stacked', 'edgecolor', 'none');
    ax = gca;
    ax.XTick = 1:length(nDataList);
    ax.XTickLabel = {'S2C 6s-TG', 'S2C 6s-AAV', 'S2C 6f-TG', 'S2C S1 6s-AAV'};
    ax.XTickLabelRotation = 30;
    ax.YTick = 0:0.1:0.7;
    box off
    ylim([0 0.7])
    xlim([0.5 length(nDataList)+0.5]);
    xlabel('PC components')
    ylabel('EV')  
    set(gca, 'TickDir', 'out')
    setPrint(8, 6, ['S2C_PCA_n_comp_' num2str(n_)], 'pdf') 
end

close all