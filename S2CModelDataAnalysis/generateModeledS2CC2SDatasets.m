addpath('../Func');
setDir;
load ([TempDatDir 'DataListS2CModel.mat']);

nData = 1; % plot raster and psth
load([TempDatDir DataSetList(nData).name '.mat'])
params  = DataSetList(nData).params;
binsize  = params.binsize;


nData = 2;
load([TempDatDir DataSetList(nData).name '.mat'])
for nCell = 3:4
    
    numYesTrial      = size(nDataSet(nData).unit_yes_trial, 1);
    numNoTrial       = size(nDataSet(nData).unit_no_trial, 1);
    numT             = size(nDataSet(nData).unit_yes_trial, 2);
    
    nUnitData        = [nDataSet(nCell).unit_yes_trial; nDataSet(nCell).unit_no_trial];
    nUnitData        = bsxfun(@minus, nUnitData, min(nUnitData, [], 2));
    nUnitData_       = nUnitData';
    
    fastData         = imagingToSpike(nUnitData_(:)');
    fastData         = reshape(fastData, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).mcmc_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).mcmc_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
    
    fr               = 1/binsize;
    spk              = mlspike(nUnitData, fr);
    t_frame                           = (0:numT)/fr;
    spks                              = spk.spk;
    fastData                          = zeros(numYesTrial+numNoTrial, numT);
    for ntrial                        = 1:length(spks)
        if ~isempty(spks{ntrial})
            fastData(ntrial, :)       = histcounts(spks{ntrial}, t_frame);
        end
    end

    nDataSet(nCell).mlspike_yes_trial = fastData(1:numYesTrial, :);
    nDataSet(nCell).mlspike_no_trial  = fastData(1+numYesTrial:end, :);
    
    [~, ~, data]     = peel_oopsi(nUnitData_(:)', fr);
    fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).peel_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).peel_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
    [~, ~, data]     = peel_nl_oopsi(nUnitData_(:)', fr);
    fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).peelnl_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).peelnl_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
    save(['S2CC2S_' DataSetList(nData).name '.mat'], 'nDataSet');
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
    
    fastData         = imagingToSpike(nUnitData_(:)');
    fastData         = reshape(fastData, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).mcmc_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).mcmc_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
    
    fr               = 1/binsize;
    spk              = mlspike(nUnitData, fr);
    t_frame                           = (0:numT)/fr;
    spks                              = spk.spk;
    fastData                          = zeros(numYesTrial+numNoTrial, numT);
    for ntrial                        = 1:length(spks)
        if ~isempty(spks{ntrial})
            fastData(ntrial, :)       = histcounts(spks{ntrial}, t_frame);
        end
    end

    nDataSet(nCell).mlspike_yes_trial = fastData(1:numYesTrial, :);
    nDataSet(nCell).mlspike_no_trial  = fastData(1+numYesTrial:end, :);
    
    [~, ~, data]     = peel_oopsi(nUnitData_(:)', fr);
    fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).peel_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).peel_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
    [~, ~, data]     = peel_nl_oopsi(nUnitData_(:)', fr);
    fastData         = reshape(data.spiketrain, numT, numYesTrial + numNoTrial);
    
    nDataSet(nCell).peelnl_yes_trial = fastData(:, 1:numYesTrial)'/binsize;    
    nDataSet(nCell).peelnl_no_trial  = fastData(:, 1+numYesTrial:end)'/binsize;
    
end

save(['S2CC2S_' fName '_.mat'], 'nDataSet');