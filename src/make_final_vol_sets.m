clear
close all

ReadPath = '../out/';
FigDir = '../fig/';
mkdir(FigDir)

load([ReadPath 'master_flux_struct.mat'],'master_flux_struct')
load([ReadPath 'height_fit_struct.mat'],'height_fit_struct')