%
% spikeToImaging.m
%
% This performs the transformation from the spiking data to the Ca++
% imaging data
% 
% Input:
% SpikingData  --- spiking time of a spike train of one cell in one trial  
%
% ----------------------------
% Output:
% CaImaging    --- Tx1 data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function CaImaging    = spikeTimeToLinearImaging(spikeTimes, timeSeriesData, params)
    
    tau_decay   = params(1);
    tau_rise    = params(2);
    intNoise    = params(3);

    Ca                        = zeros(length(spikeTimes), length(timeSeriesData));
        
    for nTrial                = 1:length(spikeTimes)
        tSpikeTimes           = spikeTimes{nTrial};
        for nSpike            = 1:length(tSpikeTimes)
            Delta_t           = max(timeSeriesData - tSpikeTimes(nSpike),0);
            Ca(nTrial,:)      = Ca(nTrial,:) + exp(-Delta_t/tau_decay).*(1-exp(-Delta_t/tau_rise));
        end
    end
    
    Ca                        = Ca + randn(size(Ca))*intNoise;    
    Ca(Ca<0)                  = 0;    
    CaImaging                 = Ca;
