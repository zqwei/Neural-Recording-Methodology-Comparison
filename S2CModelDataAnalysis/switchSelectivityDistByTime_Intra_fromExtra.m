%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load([TempDatDir 'Shuffle_Spikes.mat'])
nData = 1;
unitSpikeGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
n_cell = histcounts(unitSpikeGroup, 0:3);
spikeDataSet    = nDataSet;
sigma                         = 0.15 / DataSetList(1).params.binsize; % 200 ms
filterLength                  = 10;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse);
load ([TempDatDir 'DataListS2CModel.mat']);
m_cell = [2, 29, 4];

% cmap = cbrewer('qual', 'Set1', 3, 'cubic');

for nData      = [3 4]
    load([TempDatDir DataSetList(nData).name '.mat'])
    unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);  
    sample_cell_type = zeros(30, 3);
    for nsample = 1:30
        unitGroup_ = [];
        for ntype  = 0:2
            m_ = unitGroup(unitSpikeGroup == ntype);
            unitGroup_ = [unitGroup_; m_(randperm(n_cell(ntype+1))<=m_cell(ntype+1))];
        end
        sample_cell_type(nsample, :) = histcounts(unitGroup_, 0:3);
    end
    disp(mean(sample_cell_type)/35)
    disp(std(sample_cell_type)/35)
end

close all