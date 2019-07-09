addpath('../Func');
setDir;
file_ext = {'KL'};
xi = 0:0.001:0.5; 
color_ = [     0    0.4470    0.7410];
fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';
x_labels = {'Peakiness (KL)'};

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);
    for nData     = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat'])
        disp(DataSetList(nData).name)
        params      = DataSetList(nData).params;
        timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);
        numTimeBin          = size(nDataSet(1).unit_yes_trial, 2);
        yesProfileMatrix    = nan(length(nDataSet), numTimeBin);
        noProfileMatrix     = yesProfileMatrix;
        positivePeak        = false(length(nDataSet));
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

        timeStep  = DataSetList(nData).params.timeSeries;
        timeTag   = timePoints(2):timePoints(4)+13; % sample to response
        numTime   = length(timeTag);
        polein    = DataSetList(nData).params.polein;
        poleout   = DataSetList(nData).params.poleout;
        kl        = zeros(length(maxId), 1);
        for nNeuron    = 1:length(maxId)
            maxId_     = maxId;
            maxId_(nNeuron) = [];
            countMaxId = (hist(maxId_, 1:numTimeBin*2)+1)/size(actMat,1);
            countMaxId = countMaxId / sum(countMaxId);
            kl(nNeuron)= 1/length(countMaxId)*sum(log(countMaxId*length(countMaxId)));
        end
        for n_ = 1:1
            figure('visible', 'off');
            tmp_s = -kl(:, n_);
            fs = ksdensity(tmp_s, xi);
            f_ = max(fs);
            plot(xi, fs/f_, '-', 'color', color_, 'linewid', 2);
            set(gca,'XTick',[xi(1), xi(end)], 'YTick', [0, 1], 'TickDir', 'out')
            ylabel('Prob. density')
            xlabel(x_labels{n_})
            set(gca,'fontsize', 14)
            xlim([0, 0.5])
            ylim([0, 1])
            box off
            setPrint(8,4.5,[Result_ DataSetList(nData).name '_' file_ext{n_}], 'svg')
            close all;
        end
    end
end