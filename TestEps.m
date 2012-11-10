clc
clear all
tic
addpath('../shared_functions/')

EpsOut = ones(1,11);
EpsIn = 0:0.1:1;

%% RPGM
% InputFile = '../InputStore/RPGMLocs1Links.one'
% Options.Box = [1000 1000];
% Options.Mode = 'MaximalRange';
% Options.R = 100;
% Options.TraceMode = 'Discrete';

%% SLAW
InputFile = '/home/alireza/Applications/DownloadedProjects/SLAW/Slaw.one'
Options.InputIsContacts = false;
Options.Box = [2000 2000];
Options.Mode = 'MaximalRange';
Options.TraceMode = 'ReadAll';
Options.R = 60;
Options.Vm = 33;

%% RUN
for i = 1 : 11
    Options.OutputFile = sprintf('./Results/Eps/SLAWInferredNIEps%d.one',EpsIn(i));
    Options.epsIn = EpsIn(i);
    Options.epsOut = EpsOut(i);
    IPOptimizerWrapper( InputFile, Options );
end
toc