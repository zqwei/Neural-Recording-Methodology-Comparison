addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
params  = DataSetList(nData).params;
binsize  = params.binsize;


nData = 2;
load([TempDatDir DataSetList(nData).name '.mat'])
for nData = 3:4

    for nCell = 1:length(nDataSet)
        numYesTrial      = size(nDataSet(nCell).unit_yes_trial, 1);
        numNoTrial       = size(nDataSet(nCell).unit_no_trial, 1);
        numT             = size(nDataSet(nCell).unit_yes_trial, 2);

        nUnitData        = [nDataSet(nCell).unit_yes_trial; nDataSet(nCell).unit_no_trial];
        nUnitData        = bsxfun(@minus, nUnitData, min(nUnitData, [], 2));
        nUnitData_       = nUnitData';
        fr               = 1/binsize;

        [~, ~, data]     = peel_oopsi(nUnitData_(:)', fr);
        fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);

        nDataSet(nCell).peel_yes_trial = fastData(:, 1:numYesTrial)'/binsize;
        nDataSet(nCell).peel_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    end

    save(['S2CC2S_Peel_' DataSetList(nData).name '.mat'], 'nDataSet');
end





fName = 'FineTunedModeled_GP517';
load([TempDatDir fName '.mat'])
for nCell = 1:length(nDataSet)
    disp(nCell)
    numYesTrial      = size(nDataSet(nCell).unit_yes_trial, 1);
    numNoTrial       = size(nDataSet(nCell).unit_no_trial, 1);
    numT             = size(nDataSet(nCell).unit_yes_trial, 2);

    nUnitData        = [nDataSet(nCell).unit_yes_trial; nDataSet(nCell).unit_no_trial];
    nUnitData        = bsxfun(@minus, nUnitData, min(nUnitData, [], 2));
    nUnitData_       = nUnitData';
    fr               = 1/binsize;

    [~, ~, data]     = peel_oopsi(nUnitData_(:)', fr);
    fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);

    nDataSet(nCell).peel_yes_trial = fastData(:, 1:numYesTrial)'/binsize;
    nDataSet(nCell).peel_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;

end

save(['S2CC2S_Peel_' fName '_.mat'], 'nDataSet');                                                                                                                                                                                  
