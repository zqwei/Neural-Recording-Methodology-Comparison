% 
% Comparison based on single unit acitivity
% 
% -------------------------------------------------------------------------
% version 1.0.3
% 
% 
% 
% Check all neurons stastical properity comparing to ephys
% 
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
clc

fileList                 = {SpikeFileList; CaImagingShortDelayFastFileList; CaImagingShortDelaySlowFileList;...
                            CaImagingLongDelayFastFileList; CaImagingLongDelaySlowFileList; CaImagingShortDelaySlowVirusFileList};

nData                    = 1;
load([TempDatDir DataSetList(nData).name '.mat']);
ephysCellIndex           = [DataSetList(nData).cellinfo(:).depth]  > 110 &...
                           [DataSetList(nData).cellinfo(:).depth]  < 750;
% ephysCellDepth  = [DataSetList(nData).cellinfo(ephysCellIndex).depth];
% ephysCellML     = [DataSetList(nData).cellinfo(ephysCellIndex).ML_axis];
% [E, ephysCellIndex] = sortrows([ephysCellDepth', ephysCellML'], [2, 1]);


numUnits                  = length(nDataSet(ephysCellIndex));

% Gaussian filter for spiking data
sigma                     = 0.08 / DataSetList(nData).params.binsize; % 200 ms
filterLength              = 30;
filterStep                = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse               = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse               = filterInUse / sum (filterInUse);

pValue                    = applyFuncToCompareTrialType(nDataSet(ephysCellIndex), @pValueTTest2, filterInUse);
meanDiffValue             = applyFuncToCompareTrialType(nDataSet(ephysCellIndex), @meanDiff, filterInUse);
logPValue                 = -log(pValue);
zScores                   = sign(meanDiffValue).*logPValue;     

actMat                    = logPValue;
numT                      = size(zScores,2);

bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
bumpMat                   = actMat > bumpActThres;
bumpSize                  = ones(size(actMat,1),1);
bumpStartPoint            = nan(size(actMat,1),1);
bumpSign                  = ones(size(actMat,1),1);

for nUnit = 1: size(actMat,1)
    diffActMat                = diff([bumpMat(nUnit,:),0]);
    beginPoint                = find(diffActMat==1);
    endPoint                  = find(diffActMat==-1);
    if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
    [bumplength, bumpIndex]   = max(endPoint - beginPoint);
    if isempty(bumplength); bumplength = 0; end
    bumpSize(nUnit)           = bumplength;
    if bumplength==0
        bumpStartPoint(nUnit) = numT;
        bumpSign(nUnit)       = -1;
    else 
        bumpStartPoint(nUnit) = beginPoint(bumpIndex);
        bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
    end
end

[~, similaritySort]       = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 -1 -2]);
% similaritySort            = similaritySort(end:-1:1);
figure;
imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, zScores(similaritySort,:))
colormap(french(128,2))
caxis([-10 10])
axis xy
xlabel('Time (s)')
ylabel('Neuron Index')
box off;
% title(DataSetList(nData).name, 'interpreter', 'none')


% cThres        = 0.1;
% 
% iGroup        = [DataSetList(nData).cellinfo(ephysCellIndex).ML_axis];
% histcounts(iGroup, [unique(iGroup) inf])
% [~, ~, stats] = anova1(DataSetList(nData).params.timeSeries(bumpStartPoint), ...
%                 iGroup,'off');
% c             = multcompare(stats,'display','off');
% disp(DataSetList(nData).name)
% disp('ML dependent:')
% disp(c(c(:,end)<cThres, 1:2))
% 
% 
% iGroup        = floor(([DataSetList(nData).cellinfo(ephysCellIndex).depth]-110)/40);
% histcounts(iGroup, [unique(iGroup) inf])
% [~, ~, stats] = anova1(DataSetList(nData).params.timeSeries(bumpStartPoint), ...
%                 iGroup,'off');
% c             = multcompare(stats,'display','off');
% disp('depth dependent:')
% disp(c(c(:,end)<cThres, 1:2))
% disp('-------------------')



for nData           = 2:length(fileList)
    load([TempDatDir DataSetList(nData).name '.mat']);
    ephysCellIndex  = [DataSetList(nData).cellinfo(:).AP_axis]  > 2100 &...
                      [DataSetList(nData).cellinfo(:).AP_axis]  < 2750 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  > 1100 &...
                      [DataSetList(nData).cellinfo(:).ML_axis]  < 1900;
    numUnits                  = length(nDataSet(ephysCellIndex));

    pValue                    = applyFuncToCompareTrialType(nDataSet(ephysCellIndex), @pValueTTest2, filterInUse);
    meanDiffValue             = applyFuncToCompareTrialType(nDataSet(ephysCellIndex), @meanDiff, filterInUse);
    logPValue                 = -log(pValue);
    zScores                   = sign(meanDiffValue).*logPValue;     

    actMat                    = logPValue;
    numT                      = size(zScores,2);

    bumpActThres              = 3; % > bumpActThres considering as a bump % 3 = -log(0.05)
    bumpMat                   = actMat > bumpActThres;
    bumpSize                  = ones(size(actMat,1),1);
    bumpStartPoint            = nan(size(actMat,1),1);
    bumpSign                  = ones(size(actMat,1),1);

    for nUnit = 1: size(actMat,1)
        diffActMat                = diff([bumpMat(nUnit,:),0]);
        beginPoint                = find(diffActMat==1);
        endPoint                  = find(diffActMat==-1);
        if length(endPoint)>length(beginPoint); endPoint = endPoint(2:end); end
        [bumplength, bumpIndex]   = max(endPoint - beginPoint);
        if isempty(bumplength); bumplength = 0; end
        bumpSize(nUnit)           = bumplength;
        if bumplength==0
            bumpStartPoint(nUnit) = numT;
            bumpSign(nUnit)       = -1;
        else 
            bumpStartPoint(nUnit) = beginPoint(bumpIndex);
            bumpSign(nUnit)       = actMat(nUnit, bumpStartPoint(nUnit))/zScores(nUnit, bumpStartPoint(nUnit));
        end
    end

    [~, similaritySort]       = sortrows([bumpStartPoint, bumpSize, bumpSign], [-3 -1 -2]);
    % similaritySort            = similaritySort(end:-1:1);
    figure;
    imagesc(DataSetList(nData).params.timeSeries, 1:numUnits, zScores(similaritySort,:))
    colormap(french(128,2))
    caxis([-10 10])
    axis xy
    xlabel('Time (s)')
    ylabel('Neuron Index')
    box off;
    title(DataSetList(nData).name, 'interpreter', 'none')

    
%     iGroup        = floor(([DataSetList(nData).cellinfo(ephysCellIndex).ML_axis]-1100)/200);
%     [~, ~, stats] = anova1(DataSetList(nData).params.timeSeries(bumpStartPoint), ...
%                     iGroup,'off');
%     c             = multcompare(stats,'display','off');
%     histcounts(iGroup, [unique(iGroup) inf])
%     disp(DataSetList(nData).name)
%     disp('ML dependent:')
%     disp(c(c(:,end)<cThres, 1:2))
%     
%     disp('depth dependent:')
%     if length(unique([DataSetList(nData).cellinfo(ephysCellIndex).depth]))>1  
%         [~, ~, iGroup] = unique([DataSetList(nData).cellinfo(ephysCellIndex).depth]);
%         histcounts(iGroup, [unique(iGroup)' inf])
%         [~, ~, stats] = anova1(DataSetList(nData).params.timeSeries(bumpStartPoint), ...
%                         iGroup,'off');
%         c             = multcompare(stats,'display','off');
%         disp(c(c(:,end)<cThres, 1:2))
%     end
%     iGroup        = floor(([DataSetList(nData).cellinfo(ephysCellIndex).AP_axis]-2100)/50);
%     [~, ~, stats] = anova1(DataSetList(nData).params.timeSeries(bumpStartPoint), ...
%                     iGroup,'off');
%     c             = multcompare(stats,'display','off');
%     histcounts(iGroup, [unique(iGroup) inf])
%     disp('AP dependent:')
%     disp(c(c(:,end)<cThres, 1:2))  
%     
%     disp('-------------------')
end