function ExportToONE(X, OutputFileName, Options);
% Exports mobility matrix to ONE file
% Generic Call is ExportToOne(X, filename, Box)
%   X = The matrix to be exported
%   filename = name of the export ONE file
%   Box = size of the simulation field
    if nargin < 3
        Options = [];
    end
    N = size(X,2);
    T = size(X,3);
    if (~isfield(Options,'Box'))
        Options.Box = [1000 1000];
    end
    if (~isfield(Options,'TimeSequence'))
        Options.TimeSequence = 0:T-1;
    end
    fid = fopen(OutputFileName, 'w');
    fprintf(fid, '%.1f %.1f %.1f %.1f %.1f %.1f\n',0,T,0,...
        Options.Box(1),0,Options.Box(2));
    for t = 1:T
        for i = 1:N
            fprintf(fid, '%.1f %d %.14f %.14f\n',Options.TimeSequence(t),...
                i-1,X(1,i,t),X(2,i,t));
        end
    end
    fclose(fid);
end
