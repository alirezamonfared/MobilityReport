clc
clear all
tic
addpath('./Utilities/')


disp('This is the test file for Mobility Inference Using optimization.')
disp('To recreate the mobility traces in the paper choose one of the following traces')
disp('1 = RandomWaypoint')
disp('2 = Reference Point Group Mobility')
disp('3 = Self Similar Least Action Walks')
Response = input('Choose your synthetic trace: ')
if Response == 1
    % RWP
    InputFile = './Inputs/RWPLinks.one'
    Options.OutputFile = './Results/RWPInferred.one'
    Options.Box = [1000 1000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
    Options.R = 100;
    Options.TraceMode = 'Discrete';
elseif Response == 2
    % RPGM
    InputFile = './Inputs/RPGMLinks.one'
    Options.OutputFile = './Results/RPGMInferred.one'
    Options.Box = [1000 1000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
    Options.R = 100;
    Options.TraceMode = 'Discrete';
elseif Response == 3
    % SLAW
    InputFile = './Inputs/SLAWLocs.one'
    Options.OutputFile = './Results/SLAWInferred.one'
    Options.InputIsContacts = false;
    Options.Box = [2000 2000];
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 33;
    Options.R = 60;
    Options.TraceMode = 'ReadAll';
else
    disp('Only valid Options are 1,2, and 3')
    disp('1 = RandomWaypoint')
    disp('2 = Reference Point Group Mobility')
    disp('3 = Self Similar Least Action Walks')
end

%% MAIN RUN
Options.Mode = 'MaximalRange';

IPOptimizerWrapper( InputFile, Options );
toc