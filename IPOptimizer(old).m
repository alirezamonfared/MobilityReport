function [ X info ] = IPOptimizer( CG, X0, Options )
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
IpOptions.ipopt.hessian_approximation = 'limited-memory';
IpOptions.ipopt.mu_strategy           = 'adaptive';
%IpOptions.ipopt.tol                   = 1e-7;
IpOptions.constr_viol_tol             = 1;
IpOptions.ipopt.tol                   = Options.Vm * Options.DeltaT;
%IpOptions.max_iter                    = 300;
IpOptions.ipopt.linear_solver         = 'ma27'; % We use ma77 linear solver

IpOptions.auxdata = { X0 };
%% Setting Upper and Lower Bounds
IpOptions.lb(1:2*N) = 0;
IpOptions.ub(1:N) = Options.Box(1)*1.1; % X values == X(1:N)
IpOptions.ub(N+1:2*N) = Options.Box(2)*1.1; % Y balues == X(N+1:2N)

%% Setting Constraint Bounds
IpOptions.cl = zeros(1,N*(N+1)/2);    % N(N-1)/2 range constraints + N speed constraints
IpOptions.cu = ones(1,N*(N+1)/2)*Inf; % N(N-1)/2 range constraints + N speed constraints
% Setting constraints bounds
k = 1;
for i = 1 : N
    for j = i+1 : N
        if CG(i,j) == 1
            % packet delivery ration is between Th and 1
            sigma = Options.Threshold + (1-Options.Threshold).*rand(1);
            D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
            IpOptions.cu(k) = D;             % Lower bounds on constraints.
            k = k + 1;
        elseif CG(i,j) == 0
            % packet delivery ration is between 0 and Th
            sigma = Options.Threshold.*rand(1);
            D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
            IpOptions.cl(k) = D;             % Upper bounds on constraints.
            k = k + 1;
        else
            error('Invalid value in the connectivity graph adjacency...');
        end
    end
end

%Speed bounds
IpOptions.cl(k:k+N-1) = 0;
IpOptions.cu(k:k+N-1) = Options.Vm * Options.DeltaT;
%% 
% The callback functions.
funcs.objective         = @objective;
funcs.constraints       = @constraints;
funcs.gradient          = @gradient;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobianstructure;


%% Run IPOPT.
[ X info ] = ipopt(X0,funcs,IpOptions);
%% Convert X back to a 2xN matrix
X = [X(1:N);X(N+1:2*N)];

%% Objective
function f = objective (X, auxdata)
X0 = deal(auxdata{:});
N = 0.5*size(X0,2);
f = 0;
% for i = 1 : N
%     f = f + sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2);
% end


%% Gradient
function g = gradient (X, auxdata)
X0 = deal(auxdata{:});
N = 0.5*size(X0,2);
g = zeros(2*N,1);
% for i = 1 : N
%     g(i)   = (X(i)-X0(i))/(sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2));
%     g(i+N) = (X(i+N)-X0(i+N))/(sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2));
% end
  

%% Constraints
function c = constraints (X, auxdata)
X0 = deal(auxdata{:});
N = 0.5*size(X0,2);
k = 1;
c = zeros(N*(N+1)/2,1);
for i = 1 : N
    for j = i+1 : N
        c(k) = sqrt((X(i)-X(j))^2 + (X(i+N)-X(j+N))^2);
        k = k + 1;
    end
end
for i = 1:N
    c(k) = sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2);
    k = k + 1;
end
% disp('Total Number of constraints =')
% disp(k)

%% Jacobian
function J = jacobian (X, auxdata)
% There are N(N+1)/2 constraints and 2*N variables
% Number of nonzero enteries is:
% 4N(N-1)/2 i.e. four jacobians for each range constraint +
% 2N i.e. 2 jacobians for each speed constraint
X0 = deal(auxdata{:});
N = 0.5*size(X0,2);
%J = sparse(N*(N+1)/2, 2*N);
J = spalloc(N*(N+1)/2, 2*N, 2*N^2);
k = 1;
for i = 1 : N
    for j = i+1 : N
        J(k,i)   = (X(i)-X(j))/(sqrt((X(i)-X(j))^2 + (X(i+N)-X(j+N))^2));
        J(k,j)   = (X(j)-X(i))/(sqrt((X(i)-X(j))^2 + (X(i+N)-X(j+N))^2));
        J(k,i+N) = (X(i+N)-X(j+N))/(sqrt((X(i)-X(j))^2 + (X(i+N)-X(j+N))^2));
        J(k,j+N) = (X(j+N)-X(i+N))/(sqrt((X(i)-X(j))^2 + (X(i+N)-X(j+N))^2));
        k = k + 1;
    end
end
for i = 1:N
    J(k,i)   = (X(i)-X0(i))/(sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2));
    J(k,i+N) = (X(i+N)-X0(i+N))/(sqrt((X(i)-X0(i))^2 + (X(i+N)-X0(i+N))^2));
    k = k + 1;
end

function J = jacobianstructure (auxdata)
X0 = deal(auxdata{:});
N = 0.5*size(X0,2);
%J = sparse(N*(N+1)/2, 2*N);
J = spalloc(N*(N+1)/2, 2*N, 2*N^2);
k = 1;
for i = 1 : N
    for j = i+1 : N
        J(k,i)   = 1;
        J(k,j)   = 1;
        J(k,i+N) = 1;
        J(k,j+N) = 1;
        k = k + 1;
    end
end
for i = 1:N
    J(k,i)   = 1;
    J(k,i+N) = 1;
    k = k + 1;
end
