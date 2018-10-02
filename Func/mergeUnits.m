%
% mergeUnits.m
%
%
% ----------------------------
% Output:
%
% merged unit from unit in idx i and idx j pair
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org


function newUnit = mergeUnits(sessData, idx_i, idx_j)

    iUnit        = sessData(idx_i);
    jUnit        = sessData(idx_j);
    
    newUnit      = iUnit;
    
    numTrial     = 100;
    iYesTrial    = size(iUnit.unit_yes_trial, 1);
    jYesTrial    = size(jUnit.unit_yes_trial, 1);
    iNoTrial     = size(iUnit.unit_no_trial, 1);
    jNoTrial     = size(jUnit.unit_no_trial, 1);    
    
    for nTrial   = 1:numTrial
        newUnit.unit_yes_trial(nTrial, :) = iUnit.unit_yes_trial(ceil(rand(1)*iYesTrial), :) + jUnit.unit_yes_trial(ceil(rand(1)*jYesTrial), :);
        newUnit.unit_no_trial(nTrial, :)  = iUnit.unit_no_trial(ceil(rand(1)*iNoTrial), :) + jUnit.unit_no_trial(ceil(rand(1)*jNoTrial), :);
    end

end