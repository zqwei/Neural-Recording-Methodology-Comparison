%
% setDir.m
%
% This file sets up the basic direction information of the Ca++ imaging
% data and the spiking data
%
% version 1.0
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

warning('off', 'all')
addpath('../Func/plotFuncs')
addpath('../Func/utilFuncs')
addpath('../Func/oopsi')
addpath('../Func/MLspike/')
% set(0, 'defaultfigureVisible','off')
set(0, 'defaultaxesTickDir', 'out')
set(0, 'defaultaxesLineWidth', 1.0)

SpikingDataDir              = '../../../Data_In_Use/Dataset_Comparison/ElectrophysiologyData/';
SpikingShortNuoDir          = [SpikingDataDir 'delay1e3n/'];
SpikingShortHiDir           = [SpikingDataDir 'delay1e3h/'];
SpikingShortHiIntraDir      = [SpikingDataDir 'delay1e3h_intra/'];
SpikingLongNuoDir           = [SpikingDataDir 'delay3e0n/'];
SpikingS1NuoDir             = [SpikingDataDir 'delay1e3n_s1/'];

CaImagingDataDir            = '../../../Data_In_Use/Dataset_Comparison/ImagingData/';
CaImagingShortDelayFastDir  = [CaImagingDataDir 'delay1f4t/'];
CaImagingShortDelaySlowDir  = [CaImagingDataDir 'delay1s4t/'];
CaImagingShortDelaySlowVirusDir  = [CaImagingDataDir 'delay1s4v/'];
CaImagingLongDelayFastDir   = [CaImagingDataDir 'delay3f0t/'];
CaImagingLongDelaySlowDir   = [CaImagingDataDir 'delay3s0t/'];
CaImagingSShortDelayFastDir = [CaImagingDataDir 'delay1f2t/'];
CaImagingS1SlowVirusDir     = [CaImagingDataDir 'delay1s0v_s1/'];


% Load filenames
SpikingShortNuoFileList         = dir([SpikingShortNuoDir '*.mat']);
SpikingShortHiFileList          = dir([SpikingShortHiDir '*.mat']);
SpikingShortHiIntraFileList     = dir([SpikingShortHiIntraDir 'data*.mat']);
SpikingLongNuoFileList          = dir([SpikingLongNuoDir '*.mat']);
SpikingS1NuoFileList            = dir([SpikingS1NuoDir '*.mat']);

CaImagingShortDelayFastFileList = dir([CaImagingDataDir 'delay1f4t/*.mat']);
CaImagingShortDelaySlowFileList = dir([CaImagingDataDir 'delay1s4t/*.mat']);
CaImagingShortDelaySlowVirusFileList  = dir([CaImagingDataDir 'delay1s4v/*.mat']);
CaImagingLongDelayFastFileList  = dir([CaImagingDataDir 'delay3f0t/*.mat']);
CaImagingLongDelaySlowFileList  = dir([CaImagingDataDir 'delay3s0t/*.mat']);
CaImagingSShortDelayFastFileList= dir([CaImagingDataDir 'delay1f2t/*.mat']);
CaImagingS1SlowVirusFileList    = dir([CaImagingDataDir 'delay1s0v_s1/*.mat']);


TempDatDir                   = '../TempDat/';

if ~exist(TempDatDir, 'dir')
    mkdir(TempDatDir)
end


PlotDir                      = '../Plot/';

if ~exist(PlotDir, 'dir')
    mkdir(PlotDir)
end
