%
% obtain the fake ca++ imaging data using spikes from a list of files
%
% version 1.0
%
% Comparison list
%
% Output:
% SpikeDataSet     --- yDim x 1 cells (yDims number of neurons)
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%



function DataSetOOPSI = getFakeMLSpikeCaImagingData(spikeDataSet, fr)

    DataSetOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(spikeDataSet), 1);

    for nData          = 1:length(spikeDataSet)
        DataSetOOPSI(nData).sessionIndex     = spikeDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = spikeDataSet(nData).nUnit;
        DataSetOOPSI(nData).depth_in_um      = spikeDataSet(nData).depth_in_um;
        DataSetOOPSI(nData).AP_in_um         = spikeDataSet(nData).AP_in_um;
        DataSetOOPSI(nData).ML_in_um         = spikeDataSet(nData).ML_in_um;
        DataSetOOPSI(nData).cell_type        = spikeDataSet(nData).cell_type;

        DataSetOOPSI(nData).unit_yes_trial_index = spikeDataSet(nData).unit_yes_trial_index;
        DataSetOOPSI(nData).unit_no_trial_index  = spikeDataSet(nData).unit_no_trial_index;
        caData                                   = [spikeDataSet(nData).unit_yes_trial; spikeDataSet(nData).unit_no_trial];
        spk                                      = mlspike(caData, fr);

        numYesTrial                              = length(spikeDataSet(nData).unit_yes_trial_index);
        numNoTrial                               = length(spikeDataSet(nData).unit_no_trial_index);
        numT                                     = size(spikeDataSet(nData).unit_yes_trial, 2);
        t_frame                                  = (0:numT)/fr;
        spks                                     = spk.spk;
        fastData                                 = zeros(numYesTrial+numNoTrial, numT);
        for ntrial                               = 1:length(spks)
            if ~isempty(spks{ntrial})
                fastData(ntrial, :)              = histcounts(spks{ntrial}, t_frame);
            end
        end

        DataSetOOPSI(nData).unit_yes_trial       = fastData(1:numYesTrial, :);
        DataSetOOPSI(nData).unit_no_trial        = fastData(1+numYesTrial:end, :);
    end
