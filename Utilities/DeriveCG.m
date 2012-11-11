function  CG= DeriveCG(X, R, R2);
% Derives the connectivity graph from the locations
%   X = location matrix: d x N
%   R = trnasmission range
% Returns N x N connectivity graph
if nargin < 3
    R2 = R;
else
    R = R2;
end
    CG = squareform(pdist(X') <= R);
end