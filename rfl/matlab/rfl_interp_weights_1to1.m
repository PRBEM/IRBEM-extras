function H = rfl_interp_weights_1to1(xgrid,xhat,varargin)
% H = rfl_interp_weights_1to1(xgrid,xhat,...)
% returns H [N x Nx], sparse
% matrix of weights such that H*flux(:) interpolates the
% flux model ([Nx x 1]) into the points xhat
% xgrid is [Nx x 1] provides the x values on the x grid
% xhat is [N x 1]
% options:
% 'xbc',xbc_struct - provide boundary conditions for x interpolation
%   (see rfl_bc_init for explanation of xbc_sruct)

xbc = rfl_bc_init; % default

i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'xbc',
            i = i+1;
            xbc = rfl_bc_init(varargin{i});
        otherwise
            error('Unknown argument %s',varargin{i});
    end
    i = i+1;
end

H = rfl_lbf_eval(xgrid,xhat,xbc);
