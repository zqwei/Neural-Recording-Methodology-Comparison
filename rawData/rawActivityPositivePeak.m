%
% Peakiness
% 
% Figure 7, Wei et al., 2020
%
% 
%
% author: Ziqiang Wei
% email: weiz@janelia.hhmi.org
%
% 

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

if ~exist([PlotDir 'SingleUnitsImagescWithSort'],'dir')
    mkdir([PlotDir 'SingleUnitsImagescWithSort'])
end

for nData             = 1:length(DataSetList)
    load([TempDatDir DataSetList(nData).name '.mat'])
    plotMeanActivityImagescRasterOnlyPositivePeak(nDataSet, DataSetList(nData).params, [], []); 
    setPrint(6*2, 3*3, [PlotDir 'SingleUnitsImagescWithSort/SingleUnitsImagescRasterOnlyPositivePeak_' DataSetList(nData).name])
end

close all