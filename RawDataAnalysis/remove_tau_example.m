addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

nData   = 4;
load([TempDatDir DataSetList(nData).name '.mat'])
params  = DataSetList(nData).params;

for nNeuron = 1:4000

    trial_ = nDataSet(nNeuron).unit_no_trial';
    valid_ = sum(isnan(nDataSet(nNeuron).unit_no_trial_removeoutlier), 2)>30;
    
    if (valid_)==0
        continue
    end
    
    figure
    plot(params.timeSeries, trial_(:, ~valid_ ), '-r', 'linewid', 1)
    hold on
    plot(params.timeSeries, trial_(:, valid_ ), '-k', 'linewid', 1)

    xlim([params.timeSeries(1), params.timeSeries(end)])
    box off
    set(gca, 'tickdir', 'out')
    ylabel('DF/F')
    xlabel('Time from movements')
    setPrint(8, 6, ['tau_example_cell/tau_example_cell_' num2str(nNeuron)], 'png')
    close all
end