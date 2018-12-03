% 
% obtain synthetic spikes using ca++ imaging data from a list of files
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



function DataSetOOPSI = getSyntheticSpikeDeconvDataNLVersion(caDataSet, tau_r, tau_d, params, nlParams)
    
    inv_g              = @(p, x) p(3) - 1/p(4) * log(p(2)./(x-p(1)) - 1);
    per_cent           = 0.02;
    
    DataSetOOPSI           = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),length(caDataSet), 1);
    
    for nData          = 1:length(caDataSet)     
        DataSetOOPSI(nData).sessionIndex     = caDataSet(nData).sessionIndex;
        DataSetOOPSI(nData).nUnit            = caDataSet(nData).nUnit;
        
        yesUnitData                          = caDataSet(nData).unit_yes_trial;
        noUnitData                           = caDataSet(nData).unit_no_trial;
        minData       = min([mean(yesUnitData,1), mean(noUnitData,1)]);
        maxData       = max([mean(yesUnitData,1), mean(noUnitData,1)]);
        randCell      = ceil(rand*length(nlParams));
        paramInv      = squeeze(nlParams(randCell, :));
        paramInv(1)   = minData;
        paramInv(2)   = maxData;
        
        
        DataSetOOPSI(nData).unit_yes_trial_index = caDataSet(nData).unit_yes_trial_index;
        fastData                                 = mean(caDataSet(nData).unit_yes_trial);
        fastData                                 = bsxfun(@minus, fastData, min(fastData, [], 2));
        fastData(fastData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
        fastData(fastData < paramInv(1) + paramInv(2)*per_cent)     = paramInv(1) + paramInv(2)*per_cent; 
        fastData                                                    = inv_g(paramInv, fastData);
        fastData                                 = imagingToSpikeDeconv(fastData, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_yes_trial       = fastData;

        DataSetOOPSI(nData).unit_no_trial_index  = caDataSet(nData).unit_no_trial_index;
        fastData                                 = mean(caDataSet(nData).unit_no_trial);
        fastData                                 = bsxfun(@minus, fastData, min(fastData, [], 2));
        fastData(fastData > paramInv(1) + paramInv(2)*(1-per_cent)) = paramInv(1) + paramInv(2)*(1-per_cent); 
        fastData(fastData < paramInv(1) + paramInv(2)*per_cent)     = paramInv(1) + paramInv(2)*per_cent; 
        fastData                                                    = inv_g(paramInv, fastData);
        fastData                                 = imagingToSpikeDeconv(fastData, tau_r, tau_d, params);
        DataSetOOPSI(nData).unit_no_trial        = fastData;
    end
    
end