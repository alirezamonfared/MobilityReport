function  [CGs TimeSequence] = ImportONELinks( Filename,Options )
%IMPORTONELINKS Import a contact trace.
%   [CGs TimeSequence] = IMPORTONELINKS( Filename,Options ) Opens the ONE 
%   Link format: <Time> <CONN> <SrcID> <DstID> <UP/DOWN>, 
%   and returns CGs of dimensions N x N x tm where N is the number of the nodes
%   and tm is the maximum time. Value of each entry of CGs is 1 is the link
%   is 'UP' and 0 is the link is 'DOWN'
%       Inputs: Filename: name of the ONE file containing the link trace
%       Options: options used in the import utility
%           *NodeIDsStartFromOne: In a standard ONE format file, node IDs
%           start from 1, if the input file has node IDs starting from
%           1, set this to true
%           *TraceMode: A 'Continuous' trace reader reads the timestamps in
%           the file and reports their corresponding connectivities. It
%           uses Options.DeltaT to aggregate all contacts every DeltaT
%           seconds. A 'Discrete' trace reader, reports connectivity graphs 
%           for T=0...Tm-1 Assuming that data is updated every second. In
%           this mode, timsatmaps in the ONE file should be integers.
%           'ReadAll' is a special mode that reads all teh time instants in
%           the original file and creates one snapshot for each. Note that
%           'ReadAll' is different from 'Discrete', in the sense that it
%           does not assume that we need to report for all time instants in
%           0...T-1
%           *DiscreteTimeStartFromOne: In a standard ONE format file with
%           discrete timestamps, timestamps start from 0.0, if the input
%           file has discrete timestamps starting from 1.0, set this to
%           true. Note that this is only needed if timestamps are discrete
%           *Symmetric: If set to true, ensures a symmetric connectivity
%           graph. i.e. if i-->j then j-->i.
%           *T: the desired number of timesteps for a Continuous trace to
%           be read. The corresponding DeltaT is set as DeltaT =
%           TraceLength / T.
%           *DeltaT: the timediffenrece between cinsecutive connectivity
%           snapshpts. DeltaT is 1 for discrete mode, it can be explicity
%           set for Continuous mode, or can be implicitly set by giving the
%           number of desired timesteps. For ReadAll mode it needs not to
%           be set.
%   Warning: This utility needs the node IDs to go from 1...N or 0..N-1
%            (based on the value of Options.NodeIDsStartFromOne), if
%            node IDs are custom numbers in another range, you need to
%            pre-process your trace to make node IDs go from 0..N-1 before
%            feeding it to the utility, and remeber the mapping for future
%            purposes.
%       Output: CGs: a NxNxT connectivity graph corresponding to the
%       snapshots of the connectivity in the given contact trace.
%               TimeSequence: A 1xT vector recording the original
%               timestamps in the InputFile.

if nargin < 2
    Options = [];
end

if (~isfield(Options,'NodeIDsStartFromOne'))
    NodeOffset = 1;
else
    NodeOffset = ~Options.NodeIDsStartFromOne;
end

if (~isfield(Options,'Symmetric'))
    Symmetric = false;
else
    Symmetric = Options.Symmetric;
end

% Continuous / Dicrete / ReadAll
if (~isfield(Options,'TraceMode'))
    Options.TraceMode = 'Continuous';
end

if (~isfield(Options,'DiscreteTimeStartFromOne') && ...
        strcmp(Options.TraceMode , 'Discrete'))
    TimeOffset = 1; % Assume discrete time starts from 0
elseif (isfield(Options,'DiscreteTimeStartFromOne'))
    TimeOffset = ~Options.DiscreteTimeStartFromOne;
else
    TimeOffset = 0;
end

%% Get the last line
% T is the number of desired snapshots
% DeltaT is the time difference between consecutive snapshots

if (~isfield(Options,'DeltaT') && isfield(Options,'T'))
    % First Line
    Command = sprintf('sed q %s', Filename);
    [status,result] = system(Command);
    result = regexp(result,'\s','split');
    T1 = str2double(result(1));
    % Last Line
    Command = sprintf('sed ''$!d'' %s', Filename);
    [status,result] = system(Command);
    result = regexp(result,'\s','split');
    T2 = str2double(result(1));
    Options.DeltaT = floor((T2-T1) / Options.T);
else
    Options.DeltaT = 1;
end


%%
%Number of nodes
if (~isfield(Options,'N'))
    Options.N = 1;
end

fid = fopen(Filename);

    function Splitted = SplitString(Str)
        Parts = regexp(Str,'\s','split');
        Time = str2double(Parts{1}) + TimeOffset;
        NodeID1 = str2double(Parts{3}) + NodeOffset;
        NodeID2 = str2double(Parts{4}) + NodeOffset;
        Status = Parts{5};
        if (strcmp(Status,'UP') == 1 || strcmp(Status,'up') == 1)
            Status = 1;
        elseif (strcmp(Status,'DOWN') == 1 || strcmp(Status,'down') == 1)
            Status = 0;
        else
            Status = -1;
        end
        Splitted = [Time NodeID1 NodeID2 Status];
    end

CGs = zeros(Options.N, Options.N, 1);
TimeSequence = zeros(1);



% Continuous Mode
% In this mode, we report aggregate connectivity graphs every DeltaT
% seconds. Timestamps correspond to these moments.
if (strcmp(Options.TraceMode,'Continuous'))
    PreviousTime = -1;
    FirstLine = true;
    TimeStep = 1; % Current Timestep
    CurrentCG = zeros(Options.N,Options.N);
    while 1
        tline = fgets(fid);
        if(~ischar(tline))
            break
        end
        S = SplitString(tline);
        if (FirstLine)
            PreviousTime = S(1);
            TimeSequence(1) = PreviousTime;
            FirstLine = false;
        else
            CurrentTime = S(1);
            if (CurrentTime - PreviousTime >= Options.DeltaT)
                CGs(:,:,TimeStep) = CurrentCG;
                TimeStep = TimeStep + 1;
                TimeSequence(TimeStep) = CurrentTime;
                PreviousTime = CurrentTime;
                %CGs(:,:,TimeStep) = CGs(:,:,TimeStep-1); % Use the previous connectivity as base
            end
        end
        CurrentCG(S(2),S(3)) = S(4);
        if Symmetric
            CurrentCG(S(3),S(2)) = S(4);
        end
    end
    TimeSequence(TimeStep) = CurrentTime;
    CGs(:,:,TimeStep) = CurrentCG;
    
% Discrete Mode
% It this mode we assume that time goes from 1:T and we generate a
% Connectivity graph for every value of timestamp in this range
elseif (strcmp(Options.TraceMode,'Discrete'))
    TimeStep = 1;
    while 1
        tline = fgets(fid);
        if(~ischar(tline))
            break
        end
        S = SplitString(tline);
        TimeRead = S(1);
        while TimeRead > TimeStep
            TimeStep = TimeStep + 1;
            CGs(:,:,TimeStep) = CGs(:,:,TimeStep-1);
        end
        %TimeStep == TimeRead;
        CGs(S(2),S(3),TimeRead) = S(4);
        if Symmetric
            CGs(S(3),S(2),TimeRead) = S(4);
        end
    end
    TimeSequence = 0:S(1)-1;
    % ReadAll Mode
% In this mode, we assume a continuous trace but read all the moments in
% the trace. It is equivalent to setting Options.DeltaT to time difference
% of every two consecutive mmoments in the trace
% It is useful for discrete traces that have moments with custom timestamps
% instead of 1:T style seen in the discrete mode
elseif (strcmp(Options.TraceMode,'ReadAll'))
    TimeStep = 1;
    FirstLine = true;
    while 1
        tline = fgets(fid);
        if(~ischar(tline))
            break
        end
        S = SplitString(tline);
        if (FirstLine)
            PreviousTimeRead = S(1);
            FirstLine = false;
        end
        TimeRead = S(1);
        if TimeRead > PreviousTimeRead
            TimeSequence(TimeStep) = PreviousTimeRead;
            PreviousTimeRead = TimeRead;
            TimeStep = TimeStep + 1;
            CGs(:,:,TimeStep) = CGs(:,:,TimeStep-1);
        end
        CGs(S(2),S(3),TimeStep) = S(4);
        if Symmetric
            CGs(S(3),S(2),TimeStep) = S(4);
        end
    end
    TimeSequence(TimeStep) = TimeRead;
else
    error ('Set Options.TraceMode to either Continuous or Discrete')
end

fclose(fid);
%Options;
end



