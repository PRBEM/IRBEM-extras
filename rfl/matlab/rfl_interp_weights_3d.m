function H = rfl_interp_weights_3d(xgrid,ygrid,zgrid,xhat,yhat,zhat,varargin)
% H = rfl_interp_weights_3d(xgrid,ygrid,zgrid,xhat,yhat,zhat,...)
% returns H [N x (Nx * Ny * Nz)], sparse
% matrix of weights such that H*flux(:) interpolates the
% flux model ([Nx x Ny x Nz]) into the points xhat, yhat, zhat
% (index ordering follows ndgrid)
% xgrid is [Nx x 1] provides the x values on the x grid
% ygrid is [Ny x 1] provides the y values on the y grid
% zgrid is [Nz x 1] provides the z values on the z grid
% xhat, yhat, and zhat are each [N x 1]
% options:
% 'xbc',xbc_struct - provide boundary conditions for x interpolation
% 'ybc',ybc_struct - provide boundary conditions for y interpolation
% 'zbc',zbc_struct - provide boundary conditions for z interpolation
%   (see rfl_bc_init for explanation of *bc_sruct)

xbc = rfl_bc_init; % default
ybc = xbc;
zbc = xbc;

i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'xbc',
            i = i+1;
            xbc = rfl_bc_init(varargin{i});
        case 'ybc',
            i = i+1;
            ybc = rfl_bc_init(varargin{i});
        case 'zbc',
            i = i+1;
            zbc = rfl_bc_init(varargin{i});
        otherwise
            error('Unknown argument %s',varargin{i});
    end
    i = i+1;
end

Nx = numel(xgrid);
Ny = numel(ygrid);
Nz = numel(zgrid);

vx = sparse(rfl_lbf_eval(xgrid,xhat,xbc));
vy = sparse(rfl_lbf_eval(ygrid,yhat,ybc));
vz = sparse(rfl_lbf_eval(zgrid,zhat,zbc));
[ix,iy,iz] = ind2sub([Nx Ny Nz],1:(Nx*Ny*Nz));
H = vx(:,ix).*vy(:,iy).*vz(:,iz);
