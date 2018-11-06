%
% imagingToSpike.m
%
% This performs the transformation from the Ca++ to the spiking data
% imaging data
% 
% Input:
% nDataSet  --- nTrial x nTime
%
% ----------------------------
% Output:
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function contData = imagingToSpikeDeconv(nDataSet, tau_r, tau_d, params)

    
    contData          = nDataSet;    
    numTrial          = size(nDataSet, 1);
    binsize           = params.binsize;
    timeSeries        = params.timeSeries;
    pTime             = binsize:binsize:length(timeSeries)*binsize;
    deconv_filter     = (1 - exp(-pTime/tau_r)).* exp(-pTime/tau_d);
    deconv_filter     = deconv_filter/norm(deconv_filter);
    C                 = gallery('circul', [deconv_filter, zeros(1, length(deconv_filter))]);
    C                 = C(1:length(deconv_filter), 1:length(deconv_filter));
    C                 = C';
    invC              = inv(C);
    
    
    for nTrial        = 1:numTrial
        xdMODWT       = contData(nTrial, :);
        xdMODWT       = xdMODWT - min(xdMODWT);
        
        % XN            = contData(nTrial, :);
        % xdMODWT       = wden(XN,'modwtsqtwolog','s','mln',4,'sym4');
        contData(nTrial, :) = invC * xdMODWT'/binsize;
    end
