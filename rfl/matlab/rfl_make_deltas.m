function d = rfl_make_deltas(grid,options)
% d = rfl_make_deltas(grid,options)
% d are the default weights for a given 1-d grid

if length(grid)==1, % assume it's a tgrid: only one time means it's actually dt
    d = grid;
end

if nargin < 2,
    options = [];
end

if isfield(options,'int_method'),
    int_method = options.int_method;
else
    int_method = 'trapz'; % default
end

d = nan(size(grid));
switch(int_method),
    case 'trapz',
        d(1) = (grid(2)-grid(1))/2;
        d(2:(end-1)) = (grid(3:end)-grid(1:(end-2)))/2;
        d(end) = (grid(end)-grid(end-1))/2;
    otherwise
        error('Unknown int_method "%s"',options.int_method);
end
