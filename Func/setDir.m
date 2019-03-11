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

TempDatDir                   = '../TempDat/';

if ~exist(TempDatDir, 'dir')
    mkdir(TempDatDir)
end

url_link = https://ndownloader.figshare.com/files/14241383
https://ndownloader.figshare.com/files/13631783

websave([TempDatDir 'Shuffle_Spikes.mat'])

PlotDir                      = '../Plot/';

if ~exist(PlotDir, 'dir')
    mkdir(PlotDir)
end
