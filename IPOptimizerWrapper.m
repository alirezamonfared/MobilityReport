function XInf = IPOptimizerWrapper( InputFile, Options )
%IPOPTIMIZERWRAPPER Optimizer for mobility inference
%   XInf = IPOPTIMIZERWRAPPER ( InputFile, Options ) infers a mobility
%   trace based on the given InputFile and specified Options.
%   Inputs: InpuutFile: the name of the file containing the ONE trace for
%           the input. Note that input should be a contact trace in ONE format.(
%           the input can also be a mobility trace in ONE format. This mode is only
%           useful in simulations and to use it Optiosn.InputIsContacts needs to
%           be set to false.
%           Options: *Ipoptions: Opotions passed to Ipopt
%                    *Vm: Maximum Speed (m/s) (Default: 10)
%                    *R: Transmission Range (m) (Default: 100)
%                    *Box: 1x2 array containing the dimensions of simulaiotn
%                    field (Default [1000 1000]).
%                    *X0: Initial location of the ndoes. If not given, we
%                    will find a suitable initial point using the
%                    Mulitdimensional Scaling mechanism.
%                    *InputIsContacts: If set to false, means that input
%                    trace is a mobiltiy trace (Defaul: false).
%                    *OutputFile: Name of the file to store the output
%                    mobiltiy trace (Default: './Results/Results.one').
%                    *epsIn: Safety margin for connection constraints. We
%                    assume that connected nodes are nearer than
%                    R(1-epsIn) (Default: 0).
%                    *epsOut: Safety margin for disconnection constraints. We
%                    assume that disconnected nodes are farther than
%                    R(1+epsOut) (Default: 0).
%                    *TraceMode: Mode to read InputFile using
%                    ImportONELinks (see its help for more info) Choices
%                    are 'Continuous', 'Discrete', and 'ReadAll' (Default:
%                    Discrete).
%                    *Map: a 2xM array were the first row contains packet
%                    delivery ratios and the second row contains distances.
%                    It is used for future development of upcoming versions
%                    of our optimizer. For this version, if you call
%                    'Map' will be automatically set.
%                    *Threshold: Splitting packet delivery ratio of
%                    connected and dsiconnected nodes. It's foreseen for
%                    the furture version. Do not change in this version.
%                    (Default: 0.5).
%                    *Mode: Monility inference mode. 'MaximalRange' is the 
%                    mode used in this version. Other modes will be
%                    implemented later (default: MaximalRange).
%                    *N: Number of nodes. (This is usually guseed from the
%                    InputFile, but some reading modes require it to be
%                    explicitly set)
%                    *T: Number of timesteps. (This is usually guseed from the
%                    InputFile, but some reading modes require it to be
%                    explicitly set. For example in the 'Continuous'
%                    reading mode of contact traces, an explicit value of T
%                    corresponds the number of timesteps that one desires
%                    the extracted contact trace to have.
%                    *TimeSequence: The sequence of timestamps read from
%                    the InputFile to which enteries of the extraced
%                    connectivity graph correspond.(This is usually guseed
%                    from the InputFile, but some reading modes require it to be
%                    explicitly set).

%       

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
    Options.DistanceFilling = 'hopCount'; % FIX ME later
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

