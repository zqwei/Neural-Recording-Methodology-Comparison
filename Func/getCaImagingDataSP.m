% 
% obtain the ca++ imaging data from a list of files
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

function CaImagingDataSet    = getCaImagingDataSP(CaImagingDir, CaImagingFileList, minNumTrialToAnalysis, paramsROI)

    CaImagingDataSet   = repmat(struct('sessionIndex',1, 'nUnit', 1, ...
                                'unit_yes_trial', 1, 'unit_no_trial', 1),1000, 1);
    tot_Unit           = 0;    
    h                  = waitbar(0,'Initializing data loads...');
    for nfile = 1:length(CaImagingFileList)
        fname               = CaImagingFileList(nfile).name;
        load([CaImagingDir fname])
        numTrial        = length(s.trialTypeMat(1,:));
        unit_no_trial   = s.trialTypeMat(1,:);
        unit_yes_trial  = s.trialTypeMat(2,:);
        error_no_trial  = s.trialTypeMat(3,:);
        error_yes_trial = s.trialTypeMat(4,:);
        cue_time        = s.eventSeriesArrayHash.value{7}.eventTimes;
        pole_time       = s.eventSeriesArrayHash.value{1}.eventTimes;
        cue_trial       = s.eventSeriesArrayHash.value{7}.eventTrials;
        pole_trial      = s.eventSeriesArrayHash.value{1}.eventTrials;
        delay_dur       = zeros(numTrial, 1);
        pole_dur        = zeros(numTrial, 1);
        cue_time_       = zeros(numTrial, 1);
        for ntrial      = 1:numTrial
            if sum(pole_trial==ntrial)==2 && sum(cue_trial==ntrial)==2
                cue_    = min(cue_time(cue_trial==ntrial));
                pole_   = max(pole_time(pole_trial==ntrial));
                pole_m  = min(pole_time(pole_trial==ntrial));
                delay_dur(ntrial) = cue_ - pole_;
                pole_dur(ntrial)  = pole_ - pole_m;
                cue_time_(ntrial) = cue_;
            end
        end        
        valid_dur       = delay_dur>1040 & delay_dur<1060 & pole_dur>1010 & pole_dur<1030;
        unit_no_trial(~valid_dur)   = false;
        unit_yes_trial(~valid_dur)  = false;
        error_no_trial(~valid_dur)  = false;
        error_yes_trial(~valid_dur) = false;
        unit_no_trial_index         = find(unit_no_trial);
        unit_yes_trial_index        = find(unit_yes_trial);
        error_no_trial_index        = find(error_no_trial);
        error_yes_trial_index       = find(error_yes_trial);
        
        % dFF time series arrays
        dff_list        = s.timeSeriesArrayHash.value(2:end);
        for ndff        = 1:length(dff_list)
            dff_        = dff_list{ndff};
            ids         = dff_.ids;
            time_       = dff_.time;
            numUnit     = length(ids);
            trials      = unique(dff_.trial);
            trials      = trials(trials<=numTrial & trials>0);
%             invalid_dff = sum(isnan(dff_.valueMatrix))>0;
%             invalid_trials = unique(dff_.trial(invalid_dff));
%             trials      = setdiff(trials, invalid_trials);
            trials_index= unique(trials);
            numYesTrial = sum(unit_yes_trial(trials));
            numNoTrial  = sum(unit_no_trial(trials));
            
            if numYesTrial > minNumTrialToAnalysis && numNoTrial > minNumTrialToAnalysis
                % yes
                trial_ = intersect(trials_index,unit_yes_trial_index);
                unitTimeTrial  = zeros(numUnit, length(paramsROI.timeWindowIndexRange), length(trial_));
                for mtrial     = 1:length(trial_)
                    ntrial     = trial_(mtrial);
                    time_indx  = find(dff_.trial == ntrial);
                    time_trial = time_(time_indx);
                    [~, indx]  = min(abs(time_trial - cue_time_(ntrial)));
                    slice_     = indx + paramsROI.timeWindowIndexRange;
                    unitTimeTrial(:, :, mtrial) = dff_.valueMatrix(:, time_indx(slice_));
                end
                
                unit_yes_TimeTrial = unitTimeTrial;
                
                % no
                trial_ = intersect(trials_index,unit_no_trial_index);
                unitTimeTrial  = zeros(numUnit, length(paramsROI.timeWindowIndexRange), length(trial_));
                for mtrial     = 1:length(trial_)
                    ntrial     = trial_(mtrial);
                    time_indx  = find(dff_.trial == ntrial);
                    time_trial = time_(time_indx);
                    [~, indx]  = min(abs(time_trial - cue_time_(ntrial)));
                    slice_     = indx + paramsROI.timeWindowIndexRange;
                    unitTimeTrial(:, :, mtrial) = dff_.valueMatrix(:, time_indx(slice_));
                end
                
                unit_no_TimeTrial = unitTimeTrial;
                
                % yes error
                trial_ = intersect(trials_index,error_yes_trial_index);
                unitTimeTrial  = zeros(numUnit, length(paramsROI.timeWindowIndexRange), length(trial_));
                for mtrial     = 1:length(trial_)
                    ntrial     = trial_(mtrial);
                    time_indx  = find(dff_.trial == ntrial);
                    time_trial = time_(time_indx);
                    [~, indx]  = min(abs(time_trial - cue_time_(ntrial)));
                    slice_     = indx + paramsROI.timeWindowIndexRange;
                    unitTimeTrial(:, :, mtrial) = dff_.valueMatrix(:, time_indx(slice_));
                end
                
                error_yes_TimeTrial = unitTimeTrial;
                
                
                % no error
                trial_ = intersect(trials_index,error_no_trial_index);
                unitTimeTrial  = zeros(numUnit, length(paramsROI.timeWindowIndexRange), length(trial_));
                for mtrial     = 1:length(trial_)
                    ntrial     = trial_(mtrial);
                    time_indx  = find(dff_.trial == ntrial);
                    time_trial = time_(time_indx);
                    [~, indx]  = min(abs(time_trial - cue_time_(ntrial)));
                    slice_     = indx + paramsROI.timeWindowIndexRange;
                    unitTimeTrial(:, :, mtrial) = dff_.valueMatrix(:, time_indx(slice_));
                end
                
                error_no_TimeTrial = unitTimeTrial;
                
                for nUnit       = 1: numUnit
                    
                    yes_valid   = squeeze(~sum(isnan(unit_yes_TimeTrial(nUnit,:,:)), 2)>0);
                    no_valid    = squeeze(~sum(isnan(unit_no_TimeTrial(nUnit,:,:)), 2)>0);
                    
                    if sum(yes_valid) < minNumTrialToAnalysis || sum(no_valid) < minNumTrialToAnalysis
                        continue;
                    end
                    
                    tot_Unit    = tot_Unit + 1;
                    CaImagingDataSet(tot_Unit).sessionIndex         = nfile;
                    CaImagingDataSet(tot_Unit).nUnit                = ids(nUnit);
                    CaImagingDataSet(tot_Unit).unit_yes_trial       = squeeze(unit_yes_TimeTrial(nUnit,:,yes_valid))';
                    CaImagingDataSet(tot_Unit).unit_yes_trial_index = intersect(trials_index,unit_yes_trial_index);
                    CaImagingDataSet(tot_Unit).unit_no_trial        = squeeze(unit_no_TimeTrial(nUnit,:,no_valid))';
                    CaImagingDataSet(tot_Unit).unit_no_trial_index  = intersect(trials_index,unit_no_trial_index);

                    CaImagingDataSet(tot_Unit).unit_yes_error       = squeeze(error_yes_TimeTrial(nUnit,:,:))';
                    CaImagingDataSet(tot_Unit).unit_yes_error_index = intersect(trials_index,error_yes_trial_index);
                    CaImagingDataSet(tot_Unit).unit_no_error        = squeeze(error_no_TimeTrial(nUnit,:,:))';   
                    CaImagingDataSet(tot_Unit).unit_no_error_index  = intersect(trials_index,error_no_trial_index);
                    
                    CaImagingDataSet(tot_Unit).depth_in_um          = nan;
                    CaImagingDataSet(tot_Unit).ML_in_um             = nan;
                    CaImagingDataSet(tot_Unit).AP_in_um             = nan;
                    CaImagingDataSet(tot_Unit).cell_type            = 'non_classified';
                end
                
            end
        end
    
        waitbar(nfile/length(CaImagingFileList), h, sprintf('%d of %d files have been finished...',nfile, length(CaImagingFileList)));
    end
    
    if tot_Unit < 1000
        CaImagingDataSet = CaImagingDataSet(1:tot_Unit);
    end
    
    close (h)