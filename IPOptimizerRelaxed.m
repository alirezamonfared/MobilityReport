function [ X info ] = IPOptimizerRelaxed( CG, X0, Options )
%IPOPTIMIZER uses Ipopt package to solve one instant of the optimization
%for  a given connectivity graph
%   Detailed explanation goes here
%   Needs Options.Map, Options.Threshold, Options.Box 

% FIXME THIS IS NOT PORTABLE
% you need to change this line to the correct address where you have
% installed IpOpt. This can be also found in the file startup.m in your 
% '<DIR>/Ipopt/contrib/MatlabInterface/examples' directory
addpath /home/alireza/Applications/IpOpt/install/lib

%% Linearizing  the dxN matrix, X0
X0 = [X0(1,:) X0(2,:)];
%Extracting number of nodes
N = size(CG,1);
%% Set the IPOPT options.
if ~isfield(Options,'IpOptions')
    IpOptions.ipopt.hessian_approximation = 'limited-memory';
    IpOptions.ipopt.mu_strategy           = 'adaptive';
    %IpOptions.ipopt.tol                   = 1e-7;
    %IpOptions.ipopt.constr_viol_tol       = Options.Vm * Options.DeltaT;
    %IpOptions.ipopt.tol                   = Options.Vm * Options.DeltaT;
    IpOptions.ipopt.max_iter              = 300;
    IpOptions.ipopt.linear_solver         = 'ma27'; % We use ma77 linear solver
else
    IpOptions = Options.IpOptions;
end

%% Manually disable warm-starting 
IpOptions.ipopt.warm_start_init_point = 'no';

%% Set auxiliary data
IpOptions.auxdata = { N Options.I Options.J};
%% Setting Upper and Lower Bounds
IpOptions.lb(1:2*N) = 0;
IpOptions.ub(1:N) = Options.Box(1); % X values == X(1:N)
IpOptions.ub(N+1:2*N) = Options.Box(2); % Y balues == X(N+1:2N)

%% Setting Constraint Bounds
%IpOptions.cl = zeros(1,N*(N-1)/2);    % N(N-1)/2 range constraints + N speed constraints
%IpOptions.cl = ones(1,N*(N-1)/2)*(-Inf);    % N(N-1)/2 range constraints + N speed constraints
IpOptions.cl = ones(1,N*(N-1)/2)*(-eps);    % N(N-1)/2 range constraints + N speed constraints
IpOptions.cu = ones(1,N*(N-1)/2)*Inf; % N(N-1)/2 range constraints + N speed constraints
% Setting constraints bounds
k = 1;
Ds = Options.Map(:,2);
Sigmas = Options.Map(:,1);
for i = 1 : N
    for j = i+1 : N
        if CG(i,j) == 1
            % sigma = Options.Threshold + (1-Options.Threshold).*rand(1);
            % D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
            sigma = round((Options.Threshold + (1-Options.Threshold).*rand(1))*100)/100; % Round to two decimal places
            D = Ds(find(Sigmas == sigma, 1, 'first')); % Look up
            IpOptions.cu(k) = (D*(1-Options.epsIn)+eps)^2;             % Upper bounds on constraints.
            k = k + 1;
        elseif CG(i,j) == 0
            % packet delivery ration is between 0 and Th
            % sigma = Options.Threshold.*rand(1);
            % D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
            sigma = round((Options.Threshold.*rand(1))*100)/100;  % Round to two decimal places
            D = Ds(find(Sigmas == sigma, 1, 'first')); % Look up
            IpOptions.cl(k) = (D*(1+Options.epsOut)-eps)^2;             % Lower bounds on constraints.
            k = k + 1;
        else
            error('Invalid value in the connectivity graph adjacency...');
        end
    end
end

%% 
% The callback functions.
funcs.objective         = @objective;
funcs.constraints       = @constraints;
funcs.gradient          = @gradient;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobianstructure;
% funcs.hessian           = @hessian;
% funcs.hessianstructure  = @hessianstructure;

%% Run IPOPT.
[ X info ] = ipopt(X0,funcs,IpOptions);
%% Convert X back to a 2xN matrix
X = [X(1:N);X(N+1:2*N)];

%% Objective
function f = objective (X, auxdata)
f = 0;


%% Gradient
function g = gradient (X, auxdata)
N = auxdata{1};
g = zeros(2*N,1);
  

%% Constraints
function c = constraints (X, auxdata)
[N I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
c = zeros(N*(N-1)/2,1);
c(1:N1) = (X(I)-X(J)).^2 + (X(I+N)-X(J+N)).^2;

%% Jacobian
function J = jacobian (X, auxdata)
% There are N(N+1)/2 constraints and 2*N variables
% Number of nonzero enteries is:
% 4N(N-1)/2 i.e. four jacobians for each range constraint +
% 2N i.e. 2 jacobians for each speed constraint
[N I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
Rows = zeros(2*N^2-2*N, 1);
Cols = zeros(2*N^2-2*N, 1);
Vals = zeros(2*N^2-2*N, 1);

Rows(1:N1) = 1:N1;
Cols(1:N1) = I;
Vals(1:N1) = 2.*(X(I)-X(J));

Rows(N1+1:2*N1) = 1:N1;
Cols(N1+1:2*N1) = J;
Vals(N1+1:2*N1) = 2.*(X(J)-X(I));

Rows(2*N1+1:3*N1) = 1:N1;
Cols(2*N1+1:3*N1) = I+N;
Vals(2*N1+1:3*N1) = 2.*(X(I+N)-X(J+N));

Rows(3*N1+1:4*N1) = 1:N1;
Cols(3*N1+1:4*N1) = J+N;
Vals(3*N1+1:4*N1) = 2.*(X(J+N)-X(I+N));

J = sparse(Rows, Cols, Vals,N*(N-1)/2, 2*N, 2*N^2-2*N);

function J = jacobianstructure (auxdata)
[N I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
Rows = zeros(2*N^2-2*N, 1);
Cols = zeros(2*N^2-2*N, 1);
Vals = zeros(2*N^2-2*N, 1);

Rows(1:N1) = 1:N1;
Cols(1:N1) = I;
Vals(1:N1) = 1;

Rows(N1+1:2*N1) = 1:N1;
Cols(N1+1:2*N1) = J;
Vals(N1+1:2*N1) = 1;

Rows(2*N1+1:3*N1) = 1:N1;
Cols(2*N1+1:3*N1) = I+N;
Vals(2*N1+1:3*N1) = 1;

Rows(3*N1+1:4*N1) = 1:N1;
Cols(3*N1+1:4*N1) = J+N;
Vals(3*N1+1:4*N1) = 1;
J = sparse(Rows, Cols, Vals,N*(N-1)/2, 2*N, 2*N^2-2*N);

%% Hessian
%% Hessian
function H = hessian (X, sigma, lambda, auxdata)  
[N I J]= deal(auxdata{:});
H = sigma * zeros(2*N,2*N);
k = 1;
for i = 1 : N
    for j = i+1 : N
        Hk = zeros(2*N,2*N);
        Hk(i,i)     = 2;
        %Hk(i,j)     = 2;
        %Hk(i,i+N)   = 2;
        %Hk(i,j+N)   = 2;
        Hk(j,i)     = 2;
        Hk(j,j)     = 2;
        %Hk(j,i+N)   = 2;
        %Hk(j,j+N)   = 2;
        Hk(i+N,i)   = 2;
        Hk(i+N,j)   = 2;
        Hk(i+N,i+N) = 2;
        %Hk(i+N,j+N) = 2;
        Hk(j+N,i)   = 2;
        Hk(j+N,j)   = 2;
        Hk(j+N,i+N) = 2;
        Hk(j+N,j+N) = 2;
        %Hk = tril(Hk);
        H = H + lambda(k) * Hk;
        k = k + 1;
    end
end
H = sparse(H);


function H = hessianstructure (auxdata)
[N I J]= deal(auxdata{:});
H = zeros(2*N,2*N);
k = 1;
for i = 1 : N
    for j = i+1 : N
        H(i,i)     = 1;
        %H(i,j)     = 1;
        %H(i,i+N)   = 1;
        %H(i,j+N)   = 1;
        H(j,i)     = 1;
        H(j,j)     = 1;
        %H(j,i+N)   = 1;
        %H(j,j+N)   = 1;
        H(i+N,i)   = 1;
        H(i+N,j)   = 1;
        H(i+N,i+N) = 1;
        %H(i+N,j+N) = 1;
        H(j+N,i)   = 1;
        H(j+N,j)   = 1;
        H(j+N,i+N) = 1;
        H(j+N,j+N) = 1;
        %H = tril(H);
        k = k + 1;
    end
end
H = sparse(H);
% function H = hessian (X, sigma, lambda, auxdata)  
% [N I J]= deal(auxdata{:});
% H = sigma * sparse(2*N,2*N);
% N1 = N*(N-1)/2;
% H = H + lambda(1:N1) .* (sparse(I,I,2,2*N,2*N) + sparse(J,I,2,2*N,2*N) + ...
%     sparse(J,J,2,2*N,2*N) + sparse(I+N,I,2,2*N,2*N) +...
%     sparse(I+N,J,2,2*N,2*N) + sparse(I+N,I+N,2,2*N,2*N) + ...
%     sparse(J+N,I,2,2*N,2*N) + sparse(J+N,J,2,2*N,2*N) + ...
%     sparse(J+N,I+N,2,2*N,2*N) + sparse(J+N,J+N,2,2*N,2*N));
% H = sparse(H);
% 
% 
% function H = hessianstructure (auxdata)
% [N I J]= deal(auxdata{:});
% H = sigma * sparse(2*N,2*N);
% N1 = N*(N-1)/2;
% H = H + lambda(1:N1) .* (sparse(I,I,1,2*N,2*N) + sparse(J,I,1,2*N,2*N) + ...
%     sparse(J,J,1,2*N,2*N) + sparse(I+N,I,1,2*N,2*N) +...
%     sparse(I+N,J,1,2*N,2*N) + sparse(I+N,I+N,1,2*N,2*N) + ...
%     sparse(J+N,I,1,2*N,2*N) + sparse(J+N,J,1,2*N,2*N) + ...
%     sparse(J+N,I+N,1,2*N,2*N) + sparse(J+N,J+N,1,2*N,2*N));
% H = sparse(H);