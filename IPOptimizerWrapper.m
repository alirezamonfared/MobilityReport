function XInf = IPOptimizerWrapper( InputFile, Options )
%IPOPTIMIZERWRAPPER Accepts a contact tracec in One Format and creates a
% mobility trace as output, also in ONE format.
%   Detailed explanation goes here

%% Initializations
if (~isfield(Options, 'Box'))
    Options.Box = [1000 1000];
end

% Maximum Speed
if (~isfield(Options,'Vm'))
    Options.Vm = 10;
end

% safety margin for radio ranges
if (~isfield(Options,'epsIn'))
    Options.epsIn = 0;
end
if (~isfield(Options,'epsOut'))
    Options.epsOut = 0;
end

if (~isfield(Options,'TraceMode'))
    Options.TraceMode = 'Discrete';
    % Other option is Continuous
end

if (~isfield(Options,'Threshold'))
    Options.Threshold = 0.5;
end

if (~isfield(Options,'InputIsContacts'))
    Options.InputIsContacts = true;
end

if (~isfield(Options,'OutputFile'))
    Options.OutputFile = './Results/Results.one';
end
%% Creating a Map
if ~isfield(Options,'Mode')
    Options.Mode = 'MaximalRange';
end

if strcmp(Options.Mode,'MaximalRange')
    if ~isfield(Options,'R')
        Options.R = 100;
    end
    Options.Map(:,1) = round(linspace(0, 1, 101)*100)/100;
    Options.Map(:,2) = round(ones(1,101)*Options.R*100)/100;
elseif strcmp(Options.Mode,'UseMap')
    disp('Read From external map will be implemented soon')
end
Options = rmfield(Options,'Mode');


%% Set the IPOPT options.
MaxR = sqrt(Options.Box(1)^2+Options.Box(2)^2);
tols = [1e-6;MaxR;MaxR];

Options.IpOptions.ipopt.print_level = 0;
Options.IpOptions.ipopt.hessian_approximation = 'limited-memory';
Options.IpOptions.ipopt.mu_strategy           = 'adaptive';
Options.IpOptions.ipopt.tol                   = sum(tols);
Options.IpOptions.ipopt.constr_viol_tol       = tols(1);
Options.IpOptions.ipopt.compl_inf_tol          = tols(2);
Options.IpOptions.ipopt.dual_inf_tol           = tols(3);
Options.IpOptions.ipopt.hessian_constant = 'yes';
  Options.IpOptions.ipopt.linear_solver         = 'ma57'; % We use ma77 linear solver
    
%% Read ONE contact trace into a NxNxT connectivity graph, or from Options.CGs 
% which is a NxNxT matrix
if ~isempty(InputFile)
    if Options.InputIsContacts
        ImportParameters.TraceMode = Options.TraceMode;
        ImportParameters.Symmetric = true;
        if (isfield(Options,'N'))
            ImportParameters.N = Options.N;
        end
        if (isfield(Options,'T'))
            ImportParameters.T = Options.T;
        end
        [CGs TimeSequence] = ImportONELinks(InputFile, ImportParameters);
        Options = rmfield(Options,'TraceMode');
        
        % Extract Dimensions
        N = size(CGs,1);
        T = size(CGs,3);
    else
        [XOrg TimeSequence Nodes Box]= ImportONE(InputFile);
        T = size(XOrg,3);
        N = size(XOrg,2);
        % Create CGs
        CGs = zeros(N,N,T);
        for t = 1 : T
            CGs(:,:,t) = DeriveCG(XOrg(:,:,t),1,Options.R);
        end
    end
else
    CGs = Options.CGs;
    N = size(CGs,1);
    T = size(CGs,3);
    if ~isfield(Options,'TimeSequence')
        TimeSequence = 1 : T;
    else
        TimeSequence = Options.TimeSequence;
        Options = rmfield(Options,'TimeSequence');
    end
    Options = rmfield(Options,'CGs');
end
% Pre-allocate locations
XInf = zeros(2,N,T);

%% Find X0 using MDS-MAP, or from a given initial point
if ~isfield(Options,'X0')
    Options.DistanceFilling = 'hopCount'; % CHANGE ME FIX ME later
    X0 = MDSMAP(CGs(:,:,1), Options);
    XInf(:,:,1) = X0;
    Options = rmfield(Options,'DistanceFilling');
else
    X0 = Options.X0;
    Options = rmfield(Options,'X0');
end

%% Set Indices
[Options.I, Options.J] = GetIndices(N);

%% Find Other X s using Ipopt
[ XInf(:,:,1) info ] = IPOptimizerRelaxed( CGs(:,:,1), X0, Options );
str = sprintf('t = %d', 1);
disp(str)
disp(info)
X0 = XInf(:,:,1);
for t = 2 : T
    str = sprintf('t = %d', t);
    disp(str)
    % For each step, X0 is the solution from previous step.
    Options.DeltaT = TimeSequence(t) - TimeSequence(t-1);
    
    % warm-start / See Ipopt Options for more info.
    Options.IpOptions.ipopt.warm_start_init_point = 'yes';
    if (strcmp(Options.IpOptions.ipopt.warm_start_init_point,'yes') && t > 2)
        Options.IpOptions.zl = info.zl;
        Options.IpOptions.zu = info.zu;
        Options.IpOptions.lambda = info.lambda;
    else
        Options.IpOptions.ipopt.warm_start_init_point = 'no';
    end
    
    [ XInf(:,:,t) info ] = IPOptimizer( CGs(:,:,t), X0, Options );
    X0 = XInf(:,:,t);
    disp(info)
end


%% Export Results
if(exist(Options.OutputFile,'file'))
    delete(Options.OutputFile)
end
ExportOptions.Box = Options.Box;
ExportOptions.TimeSequence = TimeSequence;
ExportToONE(XInf, Options.OutputFile, ExportOptions);
end

