function X = MDSMAP( CG, Options )
%MDSMAP Computes a relative map matchin the connectivity graph given in CG
%   This is an implementation of the connectivity mode of MDS-MAP
%   This functino can give the coordinates in a BOX(1) x BOX(2) gird with one
%   corner on origin
%   Ref: http://dl.acm.org/citation.cfm?id=778439

if nargin < 2
    Options = [];
end

N = size(CG,1);


if (~isfield(Options, 'Box'))
    Options.Box = [];
end

% possibilities are metricstress, metricsstress, straain, and sammon
% see matlab help for details
if (~isfield(Options, 'MDSMode'))
    Options.MDSMode = 'metricstress';
end

% possibilities are hopCount useMap
if (~isfield(Options, 'DistanceFilling'))
    Options.DistanceFilling = 'hopCount';
end

if (~isfield(Options, 'Randomize'))
    Options.Randomize = false;
end

%% Finding Distance Matrix

if strcmp(Options.DistanceFilling,'hopCount')
    distances = graphallshortestpaths(sparse(CG));
    
    % we can set the unreachable nodes to have distances twice the maximum of
    % reachable nodes
    InfD = max(max(distances(isfinite(distances)))) * 2;
    % Note that distances cannot go out of the box
    if ~isempty(Options.Box)
        InfD = min(InfD, sqrt(Options.Box(1)^2+Options.Box(2)^2));
    end
    distances(distances==Inf) = InfD;
    distances = distances .* Options.R;
%     if Options.Randomize
%         distances = distances - ones(N).*Options.R + rand(N).*Options.R;
%         distances =  distances - diag(diag(distances));
%     end
elseif strcmp(Options.DistanceFilling,'useMap') && isfield(Options,'Map') ...
        && isfield(Options,'Threshold')
    % Map is a table with Sigma = Map(:,1) for packet delivery ratios
    % and D = Map(:,2) for distances
    distances = zeros(N);
    for i = 1 : N
        for j = i+1 : N
            if CG(i,j) == 1
                % packet delivery ration is between Th and 1
                sigma = Options.Threshold + (1-Options.Threshold).*rand(1);
                D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
                distances(i,j) = D;
                distances(j,i) = D;
            elseif CG(i,j) == 0
                % packet delivery ration is between 0 and Th
                sigma = Options.Threshold.*rand(1);
                D = interp1(Options.Map(:,1), Options.Map(:,2), sigma);
                distances(i,j) = D;
                distances(j,i) = D;
            else
                error('Invalid value in the connectivity graph adjacency...');
            end
        end
    end
end

%% Applying MDS
try
    X = mdscale(distances, 2, 'criterion',Options.MDSMode);
catch err
    display(err)
    X = mdscale(distances, 2, 'criterion','metricsstress'); 
end
%% Center around 0.5*Box if given
if ~isempty(Options.Box)
    X = X + repmat(Options.Box/2,size(X,1),1);
end
X = X';

end

