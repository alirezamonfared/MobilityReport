function  CG= DeriveCG(X, Mode, R);
%DERIVECG Derive connectivity graph of a mobility trace.
%   CG= DERIVECG(X, Mode, R) derives the connectivity graph in the form of
%   a NxN matrix from a given 2xN matrix of locations and a corresponding
%   transmission range.
%   Inputs:
%   *X : location matrix: d x N
%   *Mode: Set it to 1. It will be used to implement furture options.
%   *R : trnasmission range
%   Output: 
%   *CG: N x N connectivity graph
    CG = squareform(pdist(X') <= R);
end