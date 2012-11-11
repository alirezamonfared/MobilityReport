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
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
elseif Response == 2
    % RPGM
    InputFile = './Inputs/RPGMLinks.one'
    Options.OutputFile = './Results/RPGMInferredNIEps.one'
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 10;
elseif Response == 3
    % SLAW
    InputFile = './Inputs/SLAWLocs.one'
    Options.OutputFile = './Results/SLAWInferredNIEps.one'
    Options.InputIsContacts = false;
    Options.epsIn = 0.2;
    Options.epsOut = 1;
    Options.Vm = 33;
else
    disp('Only valid Options are 1,2, and 3')
end

%% MAIN RUN
% [X TimeSequence Nodes Box]= ImportONE('../InputStore/RWPLocs1.one');
% Options.X0 = X(:,:,1);

Options.Box = [1000 1000];
Options.Mode = 'MaximalRange';
Options.R = 100;
Options.TraceMode = 'Discrete';

IPOptimizerWrapper( InputFile, Options );
toc