function  [X TimeSequence Nodes Box] = ImportONE(Filename, Options);
%IMPORTONE Import a mobility trace trace.
%   [X TimeSequence Nodes Box] = IMPORTONE( InputFile,Options ) Imports a 
%   Mobility file in the ONE mobility format
%   <Time><NodeID><Location> into a dxNxTm matrix of mobility
%   Inputs: * Filename: name of the ONE file containing the mobility
%   information(filetype is detected based on extension)
%           *Options: *Mode: Can be set to 'NewIDs' to create IDs in the
%                     range 0...N-1 corresponding to every newly seen ID in teh given
%                     mobility trace or 'KeepIDs' to make exact correspondence in the
%                     enteries of the outputted 2xNxT matrix and the IDs seen in the
%                     mobiltiy trace.
%                     *NodeIDsStartFromOne; Is only used with the mode
%                     KeepIds and is used to indicate that IDs in the
%                     original trace are in range 1..N (instead of
%                     0...N-1).
%   Outputs:
%       *X: dxNxTm matrix of the locations over time
%       *TimeSequence: A vector containing the moments of time in the input
%       file
%       *Nodes: A vector containing Node IDs in ther input
%       *Box: Size of the square simulation field in 1x2 vector.


%% Parameters Initializations
if nargin < 2
    Options = [];
end

% Use 'KeepIDs' and set 'NodeIDsStartFromOne' to false if Nodes are
% discretely given in the range 0...N-1 (typical of ONE traces),
% set 'NodeIDsStartFromOne' to rue, if nodes are discretely in the
% range 1...N
if ~isfield(Options,'Mode') % Modes are KeepIDs, NewIDs
    Options.Mode = 'NewIDs';
end

if ~isfield(Options,'NodeIDsStartFromOne')
    NodeOffset = 1;
else
    NodeOffset = ~Options.NodeIDsStartFromOne;
end

%% Reading file and initializations

[Path ,Name ,Ext] = fileparts(Filename);

File = dlmread(Filename,' ');
Box = [File(1,4) File(1,6)];
File = File(2:end,:);
File = sortrows(File,1);


% Extracting the lengths
N = length(unique(File(:,2)));
if ~isfield(Options,'InputTimeSequence')
    tm = length(unique(File(:,1)));
else
    tm = length(Options.InputTimeSequence);
end
d = 2;
TimeSequence = zeros(1,tm);
Nodes = ones(1,N)*(-1);

% Initializing X
X = zeros(d,N,tm);

NN = size(File,1);

%% Improting
t = 1;
nr = 1;
if ~isfield(Options,'InputTimeSequence')
    PreviousTime = File(1,1);
else
    PreviousTime = Options.InputTimeSequence(1);
end
TimeSequence(1) = PreviousTime;
for i = 1:NN
    CurrentNode = File(i,2);
    CurrentTime = File(i,1);
    % see if there is a new Time to record
    if isfield(Options,'InputTimeSequence') && ...
            isempty(find(Options.InputTimeSequence==CurrentTime,1))
        continue
    else
        if CurrentTime ~= PreviousTime
            t = t + 1;
            TimeSequence(t) = CurrentTime;
            PreviousTime = CurrentTime;
        end
        % see if there is a new node to record
        if strcmp(Options.Mode, 'NewIDs')
            if ~any(CurrentNode == Nodes)
                Nodes(nr) = CurrentNode;
                n = nr;
                nr = nr + 1;
            else
                n = find(Nodes == CurrentNode);
            end
        elseif strcmp(Options.Mode, 'KeepIDs')
            n = CurrentNode + NodeOffset;
        else
            error('Invalid Mode chosen')
        end
        X(:,n,t) = File(i,3:3+d-1);
    end
end

if strcmp(Options.Mode, 'KeepIDs')
   Nodes = (0:(N-1)) + NodeOffset;
end







end
