function d = rfl_make_deltas(grid,options)
% d = rfl_make_deltas(grid,options)
% d are the default weights for a given 1-d grid

d = nan(size(grid));
d(1) = grid(2)-grid(1);
d(2:(end-1)) = (grid(3:end)-grid(1:(end-2)))/2;
d(end) = grid(end)-grid(end-1);
