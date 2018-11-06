%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw activity of all neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

ylabels                 = {'Fring Rate (Hz)', 'dF/F', 'dF/F', 'dF/F', 'dR/R' };
yAxes_set               = [0 60; -0.5 2.0; -0.5 2.0; -0.5 2.0; -0.5 2.0 ; -0.5 2.0];
lowFiringThres          = [15, 0.3, 0.3, 0.3, 0.3, 0.3];

if ~exist([PlotDir 'SingleUnitsImagescWithSort'],'dir')
    mkdir([PlotDir 'SingleUnitsImagescWithSort'])
end

for nData             = [1 3 4 10]
    if nData   == 1
            load([TempDatDir DataSetList(nData).name '.mat'])
        else
            load([TempDatDir DataSetList(nData).name '_withOLRemoval.mat'])
    end
%     positivePeak = plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(nData).params, [], []);
    plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(nData).params, [], []); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeak_' DataSetList(nData).name])
end



% close all;