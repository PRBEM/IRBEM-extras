function bc = rfl_bc_init(bc)
% bc = rfl_bc_init(bc)
% an rfl boundary condition (used by interp routines)
% bc has fields:
% bc.period = period; 
%   coordinate is periodic with specified period (0 or absent for
%   non-periodic)
% bc.left = 'zero'; % clips values below first grid point to zero
% bc.left = 'extrap'; % extrapolates linear below first grid point
%   e.g., for x(1)-dx < xhat < x(1), gives (xhat-(x(1)-dx))/dx
%    where dx = x(2)-x(1)
% bc.right takes same values as .left with corresponding meaning 
%  for last points beyond grid point
% a cell array can be supplied instead and will be converted to a struct
% pass bc = [] to get default bc (clip outside grid)

if (nargin < 1) || isempty(bc),
    bc = struct; % empty 1x1 stucture
elseif iscell(bc),
    bc = struct(bc{:});
end

if ~isstruct(bc),
    error('Could not convert bc from class %s to struct',class(bc));
end

if ~isfield(bc,'period'),
    bc.period = 0;
end

if ~isfield(bc,'left'),
    bc.left = 'zero';
end

if ~isfield(bc,'right'),
    bc.right = 'zero';
end

