 function [ X info ] = IPOptimizer( CG, X0, Options )
%IPOPTIMIZER Location from connectivity problem solver
%   [ X info ] = IPOPTIMIZER ( CG, X0, Options ) performs non-metric multidimensional scaling on the
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

IpOptions.auxdata = { N X0 Options.I Options.J};
%% Setting Upper and Lower Bounds
IpOptions.lb(1:2*N) = 0;
IpOptions.ub(1:N) = Options.Box(1); % X values == X(1:N)
IpOptions.ub(N+1:2*N) = Options.Box(2); % Y balues == X(N+1:2N)

%% Setting Constraint Bounds
%IpOptions.cl = zeros(1,N*(N+1)/2);    % N(N-1)/2 range constraints + N speed constraints
%IpOptions.cl = ones(1,N*(N+1)/2)*(-Inf);    % N(N-1)/2 range constraints + N speed constraints
IpOptions.cl = ones(1,N*(N+1)/2)*(-eps);    % N(N-1)/2 range constraints + N speed constraints
IpOptions.cu = ones(1,N*(N+1)/2)*Inf; % N(N-1)/2 range constraints + N speed constraints
% Setting constraints bounds
k = 1;
Ds = Options.Map(:,2);
Sigmas = Options.Map(:,1);
for i = 1 : N
    for j = i+1 : N
        if CG(i,j) == 1
            % packet delivery ration is between Th and 1
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

%Speed bounds
IpOptions.cl(k:k+N-1) = -eps;
IpOptions.cu(k:k+N-1) = (Options.Vm * Options.DeltaT+eps)^2;
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
[N X0 I J]= deal(auxdata{:});
I2 = 1:N;
f = f + sum((X(I2)-X0(I2)).^2 + (X(I2+N)-X0(I2+N)).^2);

%% Gradient
function g = gradient (X, auxdata)
N = auxdata{1};
g = zeros(2*N,1);
[N X0 I J]= deal(auxdata{:});
I1 = 1:2*N;
g(I1) = 2.*(X(I1)-X0(I1));

%% Constraints
function c = constraints (X, auxdata)
[N X0 I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
I2 = 1:N;
c = zeros(N*(N+1)/2,1);
c(1:N1) = (X(I)-X(J)).^2 + (X(I+N)-X(J+N)).^2;
c(N1+1:N1+N) = (X(I2)-X0(I2)).^2 + (X(I2+N)-X0(I2+N)).^2;


%% Jacobian
function J = jacobian (X, auxdata)
% There are N(N+1)/2 constraints and 2*N variables
% Number of nonzero enteries is:
% 4N(N-1)/2 i.e. four jacobians for each range constraint +
% 2N i.e. 2 jacobians for each speed constraint
[N X0 I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
I2 = 1:N;
Rows = zeros(2*N^2, 1);
Cols = zeros(2*N^2, 1);
Vals = zeros(2*N^2, 1);

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

Rows(4*N1+1:4*N1+N) = N1+1:N1+N;
Cols(4*N1+1:4*N1+N) = I2;
Vals(4*N1+1:4*N1+N) = 2.*(X(I2)-X0(I2));

Rows(4*N1+N+1:4*N1+2*N) = N1+1:N1+N;
Cols(4*N1+N+1:4*N1+2*N) = I2+N;
Vals(4*N1+N+1:4*N1+2*N) = 2.*(X(I2+N)-X0(I2+N));

J = sparse(Rows, Cols, Vals,N*(N+1)/2, 2*N, 2*N^2);

function J = jacobianstructure (auxdata)
[N X0 I J]= deal(auxdata{:});
N1 = N*(N-1)/2;
I2 = 1:N;
Rows = zeros(2*N^2, 1);
Cols = zeros(2*N^2, 1);
Vals = zeros(2*N^2, 1);

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

Rows(4*N1+1:4*N1+N) = N1+1:N1+N;
Cols(4*N1+1:4*N1+N) = I2;
Vals(4*N1+1:4*N1+N) = 1;

Rows(4*N1+N+1:4*N1+2*N) = N1+1:N1+N;
Cols(4*N1+N+1:4*N1+2*N) = I2+N;
Vals(4*N1+N+1:4*N1+2*N) = 1;

J = sparse(Rows, Cols, Vals,N*(N+1)/2, 2*N, 2*N^2);

%% Hessian
function H = hessian (X, sigma, lambda, auxdata)  
[N X0 I J]= deal(auxdata{:});
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
for i = 1:N
    Hk = zeros(2*N,2*N);
    Hk(i,i)     = 2;
    %Hk(i,i+N)   = 2;
    Hk(i+N,i)   = 2;
    Hk(i+N,i+N) = 2;
    %Hk = tril(Hk);
    H = H + lambda(k) * Hk;
    k = k + 1;
end
H = sparse(H);


function H = hessianstructure (auxdata)
[N X0 I J]= deal(auxdata{:});
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
for i = 1:N
    H(i,i)     = 1;
    %H(i,i+N)   = 1;
    H(i+N,i)   = 1;
    H(i+N,i+N) = 1;
    %H = tril(H);
    k = k + 1;
end
H = sparse(H);
