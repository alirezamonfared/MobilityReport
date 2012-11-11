clc
clear all
tic
addpath('./Utilities/')

%% RWP
InputFile = './Inputs/RWPLinks.one'
Options.OutputFile = './Results/RWPInferredNIEps.one'
Options.epsIn = 0.2;
Options.epsOut = 1;
Options.Vm = 10;

%% RPGM
% InputFile = './Inputs/RPGMLinks.one'
% Options.OutputFile = './Results/RPGMInferredNIEps.one'
% Options.epsIn = 0.2;
% Options.epsOut = 1;
% Options.Vm = 10;

%% SLAW
% InputFile = './Inputs/SLAWLocs.one'
% Options.OutputFile = './Results/SLAWInferredNIEps.one'
% Options.InputIsContacts = false;
% Options.epsIn = 0.2;
% Options.epsOut = 1;
% Options.Vm = 33;

%% GaussMarkov
% InputFile = '../InputStore/GMLocs1.one'
% Options.OutputFile = './Results/GMInferredNIEps.one'
% Options.InputIsContacts = false;
% Options.epsIn = 0.2;
% Options.epsOut = 1;
% Options.Vm = 10;

%% MAIN RUN
% [X TimeSequence Nodes Box]= ImportONE('../InputStore/RWPLocs1.one');
% Options.X0 = X(:,:,1);

Options.Box = [1000 1000];
Options.Mode = 'MaximalRange';
Options.R = 100;
Options.TraceMode = 'Discrete';

IPOptimizerWrapper( InputFile, Options );
toc