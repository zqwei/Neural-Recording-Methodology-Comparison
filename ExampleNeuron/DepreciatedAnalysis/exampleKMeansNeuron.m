%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% profileAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exampleKMeansNeuron
    addpath('../Func');
    setDir;
    load ([TempDatDir 'DataListShuffle.mat']);
    % resample all data with a larger binsize
    timeBinLength = 3;
    k             = 64;
    m             = sqrt(k);
    
    if ~exist([PlotDir 'ActivityProfileKMeans'],'dir')
        mkdir([PlotDir 'ActivityProfileKMeans'])
    end
    
    for nData     = [1 3 4]
        timeStep  = DataSetList(nData).params.timeSeries(ceil(timeBinLength/2):timeBinLength:end-1);
        polein    = DataSetList(nData).params.polein;
        poleout   = DataSetList(nData).params.poleout;
        load([TempDatDir DataSetList(nData).name '.mat'])
        actMat    = computeMeanActivityMatrix(nDataSet, timeBinLength);
        [idx, ~, ~, distCenter] = kmeans(actMat, k);
        [~, cell_id] = min(distCenter);
        trialSize = size(actMat, 2)/2;
        figure;
        for nn    = 1:k
            subplot(m, m, nn)
            hold on
            plot(timeStep, actMat(cell_id(nn), 1:trialSize), '-', 'linewid', 2.0, 'color', [0 0 0.7])
            plot(timeStep, actMat(cell_id(nn), 1+trialSize:end), '-', 'linewid', 2.0, 'color', [0.7 0 0])
            xlabel('Normalized activity')
            ylabel('Time')
            xlim([timeStep(1) timeStep(end)])
            ylim([0 1])
            gridxy ([polein, poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
            title(['%' num2str(sum(idx==nn)/length(idx)*100, '%.2f') ' of pop.'])
            hold off
            box off
        end
        setPrint(8*m, 6*m, [PlotDir 'ActivityProfileKMeans/ExampleNeuronActivityProfileKMeans_K_' num2str(k, '%03d') '_' DataSetList(nData).name])
    end
    
    close all
end


function actMat      = computeMeanActivityMatrix(nDataSet, timeBinLength)
    numTimeBin       = floor(size(nDataSet(1).unit_yes_trial, 2)/timeBinLength);
    trunctSize       = mod(size(nDataSet(1).unit_yes_trial, 2), timeBinLength);
    yesProfileMatrix = nan(length(nDataSet), numTimeBin);
    noProfileMatrix  = yesProfileMatrix;
    
    for nUnit        = 1:length(nDataSet)
        yesData      = mean(nDataSet(nUnit).unit_yes_trial(:, 1:end-trunctSize));
        noData       = mean(nDataSet(nUnit).unit_no_trial(:, 1:end-trunctSize));
        maxData      = max([yesData, noData]);
        minData      = min([yesData, noData]);
        yesData      = (yesData - minData)/(maxData - minData);
        yesData      = reshape(yesData, [timeBinLength, numTimeBin]);
        yesProfileMatrix(nUnit, :) = mean(yesData);
%         maxData      = max([yesData, noData]);
%         minData      = min([yesData, noData]);
        noData       = (noData - minData)/(maxData - minData);
        noData       = reshape(noData, [timeBinLength, numTimeBin]);
        noProfileMatrix(nUnit, :) = mean(noData);
    end
    
    actMat  = [yesProfileMatrix, noProfileMatrix];
    
end
