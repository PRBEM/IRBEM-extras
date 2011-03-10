function H = rfl_interp_weights_2d(xgrid,ygrid,xhat,yhat,varargin)
% H = rfl_interp_weights_2d(xgrid,ygrid,xhat,yhat,...)
% returns H [N x (Nx * Ny)], sparse
% matrix of weights such that H*flux(:) interpolates the
% flux model ([Nx x Ny]) into the points xhat, yhat
% xgrid is [Nx x 1] provides the x values on the x grid
% ygrid is [Ny x 1] provides the y values on the y grid
% xhat and yhat are each [N x 1]
% options:
% 'xbc',xbc_struct - provide boundary conditions for x interpolation
% 'ybc',ybc_struct - provide boundary conditions for y interpolation
%   (see rfl_bc_init for explanation of *bc_sruct)

xbc = rfl_bc_init; % default
ybc = xbc;

i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'xbc',
            i = i+1;
            xbc = rfl_bc_init(varargin{i});
        case 'ybc',
            i = i+1;
            ybc = rfl_bc_init(varargin{i});
        otherwise
            error('Unknown argument %s',varargin{i});
    end
    i = i+1;
end

Nx = numel(xgrid);
Ny = numel(ygrid);
N = numel(xhat);

H = sparse(N,Nx*Ny);
vx = rfl_lbf_eval(xgrid,xhat,xbc);
vy = rfl_lbf_eval(ygrid,yhat,ybc);
for i = 1:(Nx*Ny),
    [ix,iy,iz] = ind2sub([Nx Ny],i);
    H(:,i) = vx(:,ix).*vy(:,iy);
end

