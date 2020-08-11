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
set(0, 'defaultaxesTickDir', 'out')
set(0, 'defaultaxesLineWidth', 1.0)

TempDatDir                   = '../TempDat/';

if ~exist(TempDatDir, 'dir')
    disp('Please download data from https://ndownloader.figshare.com/articles/12786296/versions/1 and save it TempDat folder')
    mkdir(TempDatDir)
end

PlotDir                      = '../Plot/';

if ~exist(PlotDir, 'dir')
    mkdir(PlotDir)
end
