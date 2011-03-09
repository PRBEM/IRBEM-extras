function v = rfl_lbf_eval(xgrid,xhat,xbc)
% v = rfl_lbf_eval(xgrid,xhat,xbc)
% xgrid is list of grid x values [Nx x 1]
% xhat is list of query x values [N x 1]
% xbc (optional) is a boundary conditions structure (see rfl_bc_init)
% v is the value of the linear basis functions for the xgrid at points xhat
% v is [N x Nx]

if nargin < 3,
    xbc = [];
end

Nx = numel(xgrid);
if Nx < 2,
    error('xgrid must contain at least 2 points');
end
xbc = rfl_bc_init(xbc);
N = numel(xhat);

v = sparse(N,Nx);
for i = 2:Nx,
    f = ((xhat >= xgrid(i-1)) & (xhat <= xgrid(i))); % allow >= and <= to include xhat == xgrid(1)
    if any(f),
        v(f,i) = (xhat(f)-xgrid(i-1))/(xgrid(i)-xgrid(i-1));
        v(f,i-1) = 1-v(f,i); % allocate rest to point one point to the left
    end
end

if isfield(xbc,'period') && xbc.period, % wrap around
    f = ((xhat > xgrid(end)-xbc.period) & (xhat <= xgrid(1))); % before first point
    if any(f),
        v(f,1) = (xhat(f)-xgrid(end)+xbc.period)/(xgrid(1)-xgrid(end)+xbc.period);
        v(f,end) = 1-v(f,1);
    end
    f = ((xhat > xgrid(end)) & (xhat < xgrid(1)+xbc.period)); % after last point
    if any(f),
        v(f,end) = (xhat(f)-xgrid(end))/(xgrid(1)+xbc.period-xgrid(end));
        v(f,1) = 1-v(f,end);
    end
else % check for extrap bc
    if strcmp(xbc.left,'extrap'), % left extrapolate
        dx = xgrid(2)-xgrid(1);
        f = ((xhat > xgrid(1)-dx) & (xhat < xgrid(1)));
        if any(f),
            v(f,1) = (xhat(f)-xgrid(1)+dx)/dx;
        end
    end
    if strcmp(xbc.right,'extrap'), % right extrapolate
        dx = xgrid(end)-xgrid(end-1);
        f = ((xhat > xgrid(end)) & (xhat < xgrid(end)+dx));
        if any(f),
            v(f,end) = (xgrid(end)+dx-xhat(f))/dx;
        end
    end
end
