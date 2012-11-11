clc
clear all
tic
addpath('./Utilities/')


disp('This is the test file for Mobility Inference Using optimization.')
disp('To recreate the mobility traces in the paper choose one of the following traces')
dis('1 = RandomWaypoint')
dis('2 = Reference Point Group Mobility')
dis('3 = Self Similar Least Action Walks')
Response = input('Choose your synthetic trace: ')
if Response == 1
    % RWP
    InputFile = './Inputs/RWPLinks.one'
    Options.OutputFile = './Results/RWPInferredNIEps.one'
    Options.Box = [1000 1000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
    Options.R = 100;
elseif Response == 2
    % RPGM
    InputFile = './Inputs/RPGMLinks.one'
    Options.OutputFile = './Results/RPGMInferredNIEps.one'
    Options.Box = [1000 1000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
    Options.R = 100;
elseif Response == 3
    % SLAW
    InputFile = './Inputs/SLAWLocs.one'
    Options.OutputFile = './Results/SLAWInferredNIEps.one'
    Options.InputIsContacts = false;
    Options.Box = [2000 2000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 33;
    Options.R = 60;
else
    disp('Only valid Options are 1,2, and 3')
    dis('1 = RandomWaypoint')
    dis('2 = Reference Point Group Mobility')
    dis('3 = Self Similar Least Action Walks')
end

%% MAIN RUN
Options.Mode = 'MaximalRange';
Options.TraceMode = 'Discrete';

IPOptimizerWrapper( InputFile, Options );
toc