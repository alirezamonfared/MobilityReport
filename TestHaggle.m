clc
clear all
tic
addpath('../shared_functions/')


InputFile = '/home/alireza/workspace/MobilityProject/MATLAB/InputStore/Haggle3-Infocom5.csv'
Options.OutputFile = './Results/HaggleInferred.one'
Options.epsIn = 0.2;
Options.epsOut = 1;
Options.Vm = 10;

%% MAIN RUN
Options.TraceMode = 'Continuous';
Options.T = 1000;
Options.N = 41;

Options.Box = [1000 1000];
Options.Mode = 'MaximalRange';
Options.R = 100;

IPOptimizerWrapper( InputFile, Options );
toc