addpath('../Func');
setDir;
file_ext = {'Peakiness'};
xi = 0.0:0.01:2.0; 
color_ = [     0    0.4470    0.7410];
fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);

    for nData     = 1:length(DataSetList)
        load([TempDatDir DataSetList(nData).name '.mat'])
        params      = DataSetList(nData).params;
        timePoints  = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

        numTimeBin          = size(nDataSet(nData).unit_yes_trial, 2);
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

        countMaxId = hist(maxId, 1:numTimeBin*2)/size(actMat,1)*100;
        [bootstat,bootsam] = bootstrp(1000,@std,countMaxId([timeTag, timeTag+numTimeBin]));
        for n_ = 1:1
            figure
            tmp_s = bootstat(:, n_);
            fs = ksdensity(tmp_s, xi);
            f_ = max(fs);
            plot(xi, fs/f_, '-', 'color', color_, 'linewid', 2);
            set(gca,'XTick',[0, 2], 'YTick', [], 'TickDir', 'out','Ycolor','none')
            box off
            setPrint(8,3,[Result_ DataSetList(nData).name '_' file_ext{n_}], 'svg')
            close all;
        end
    end
end
