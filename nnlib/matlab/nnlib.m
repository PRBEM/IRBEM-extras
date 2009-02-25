function varargout = nnlib(what,varargin)
%
% nnlib('save_training_set',filename,X,Y,s);
%
% [X,Y,s] = nnlib('load_training_set',filename);
%
% Y = nnlib('eval',X,Nh,theta,xbar,ybar,sx,sy)
% [Y,dY] = nnlib('eval',X,Nh,theta,xbar,ybar,sx,sy,dY_flag,theta_cov)
%    dY_flag & 3 = 0: don't compute dY (theta_cov ignored)
%    dY_flag & 3 = 1: for 1 x Ny dY (standard errors)
%    dY_flag & 3 = 2: for Ny x Ny dY (covariance matrix)
%    dY_flag & 4 = 4: theta_cov is really theta_hess
% [ell,grad,hess] = nnlib('ell_debug',X,Nh,Y,s,theta,xbar,ybar,sx,sy)
%
% nnlib('save_net',filename,Nx,Ny,theta,xbar,ybar,sx,sy,theta_cov);
% nnlib('save_net',filename,Nx,Ny,theta,xbar,ybar,sx,sy,theta_hess,'hess_theta');
%
% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy,theta_cov] = nnlib('load_net',filename); 
%  (theta_cov returns empty if not in file)
% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy,theta_hess] = nnlib('load_net',filename,'theta_hess'); 
%  (theta_hess returns empty if not in file)
%  NOTE: theta_cov = inv(theta_hess)
%
% Ntheta = nnlib('Ntheta',Nx,Nh,Ny)
% Nh = nnlib('Nh',Nx,Ny,Ntheta)
%
% [theta,xbar,ybar,sx,sy,cov_theta] = nnlib('fit',X,Nh,Y,s,MaxIter,epsabs,...);
% [theta,xbar,ybar,sx,sy,cov_theta] = nnlib('fit',X,initial_theta,Y,s,MaxIter,epsabs,...);
% options: 
% 'do_cov' - compute cov_theta
% 'do_hess' - compute hess_theta, return in place of cov_theta
% 'bfgs' - use vector_bfgs minimizer (default)
% 'fr' - use conjugate_fr minimizer
% 'pr' - use conjugate_pr minimizer
% 'nm' - use neldear-mead simplex


if ~libisloaded('nnlib'),
    libfile = which('nnlib.dll');
    hfile = which('nnlib.h');
    loadlibrary(libfile,hfile,'alias','nnlib');
end

switch(lower(what)),
    case {'save_training_set'},
        varargout = save_training_set(varargin{:});
    case {'load_training_set'},
        varargout = load_training_set(varargin{:});
    case {'fit'},
        varargout = fit(varargin{:});        
    case {'ell_debug'},
        varargout = ell_debug(varargin{:});        
    case {'save_net'},
        varargout = save_net(varargin{:});
    case {'load_net'},
        varargout = load_net(varargin{:});
    case {'eval'},
        varargout = eval(varargin{:});
    case {'ntheta'},
        varargout = calc_Ntheta(varargin{:});
    case {'nh'},
        varargout = calc_Nh(varargin{:});
    otherwise
        error('%s: Unknown "what": "%s"',mfilename,what);
end

if length(varargout)>nargout,
    varargout = varargout(1:nargout);
end

function out = calc_Ntheta(Nx,Nh,Ny)
Ntheta = Ny+Ny*Nh+Nh+Nh*Nx;
out = {Ntheta};

function out = calc_Nh(Nx,Ny,Ntheta)
Nh = (Ntheta-Ny)/(Ny+Nx+1);
out = {Nh};

function out = save_training_set(varargin)
filename = varargin{1};
X = varargin{2};
Y = varargin{3};
s = varargin{4};

Nt = size(X,1);
Nx = size(X,2);
Ny = size(Y,2);

if Ny==1,
    if numel(s)==1, % scalar, 1 x 1
        flag = 0;
    elseif max(size(s))==Nt, % Nt x 1 or Nt x 1 x 1
        flag = 2;
    else % 1 x Ny
        flag = 1;
    end
else
    flag = sum(size(s)>1);
end
s = permute(s,ndims(s):-1:1);
calllib('nnlib','nnlib_save_training_set',filename,Nt,Nx,X',Ny,Y',s,flag);
out = {};

function out = load_training_set(varargin)
filename = varargin{1};
if ~exist(filename,'file'),
    error('%s: load_training_set, file not found "%s"',mfilename,filename);
end

Ntptr = libpointer('uint32Ptr',0);
Nxptr = libpointer('uint32Ptr',0);
Nyptr = libpointer('uint32Ptr',0);
flagptr = libpointer('uint32Ptr',99999);
Xptr = libpointer('doublePtr',nan);
Yptr = libpointer('doublePtr',nan);
sptr = libpointer('doublePtr',nan);
% call first time with flag=99999 to get sizes
calllib('nnlib','nnlib_load_training_set',filename,Ntptr,Nxptr,Xptr,Nyptr,Yptr,sptr,flagptr);
Nx = get(Nxptr,'value');
Ny = get(Nyptr,'value');
Nt = get(Ntptr,'value');
flag = get(flagptr,'value');
Xptr = libpointer('doublePtr',nan(Nx,Nt));
Yptr = libpointer('doublePtr',nan(Ny,Nt));
switch(bitand(flag,3)),
    case 0,
        size_s = [1 1];
    case 1,
        size_s = [Ny 1];
    case 2,
        size_s = [Ny Nt];
    case 3,
        size_s = [Ny Ny Nt];
end
sptr = libpointer('doublePtr',nan(size_s));
% call again to actually read data
calllib('nnlib','nnlib_load_training_set',filename,Ntptr,Nxptr,Xptr,Nyptr,Yptr,sptr,flagptr);

X = get(Xptr,'value')';
Y = get(Yptr,'value')';
s = reshape(get(sptr,'value'),size_s); % libpointers can only hold matrices, no 3-D stuff
s = permute(s,ndims(s):-1:1);

out = {X,Y,s};

function out = eval(varargin)
X = varargin{1};
Nt = size(X,1);
Nx = size(X,2);
Nh = varargin{2};
theta = varargin{3};
Ntheta = length(theta);
Ny = (Ntheta-(Nx+1)*Nh)/(Nh+1);
xbar = varargin{4};
ybar = varargin{5};
sx = varargin{6};
sy = varargin{7};

Yptr = libpointer('doublePtr',nan(Ny,Nt));
if length(varargin)>=8,
    dY_flag = varargin{8};
    theta_cov = varargin{9};
    if bitand(dY_flag,3) == 0,
        dYptr = libpointer('doublePtr',nan); % don't need this
    elseif bitand(dY_flag,3) == 1,
        dYptr = libpointer('doublePtr',nan(Ny,Nt));
    else % bitand(dY_flag,3) == 2,
        dYptr = libpointer('doublePtr',nan(Ny*Ny,Nt));
    end
else
    dY_flag = 0;
    theta_cov = 0;
    dYptr = libpointer('doublePtr',nan);
end

calllib('nnlib','nnlib_eval',Nt,Nx,X',Nh,theta,xbar,ybar,sx,sy,Ny,Yptr,dY_flag,theta_cov',dYptr);

Y = get(Yptr,'value')';
if bitand(dY_flag,3),
    dY = get(dYptr,'value')';
    if bitand(dY_flag,3) == 2,
        dY = reshape(dY,[Nt,Ny,Ny]);
    end
    out = {Y,dY};
else
    out = {Y};
end

function out = save_net(varargin)
% nnlib('save_net',filename,Nx,Ny,theta,theta_cov,hess_theta_flag);
filename = varargin{1};
Nx = varargin{2};
Ny = varargin{3};
theta = varargin{4};
xbar = varargin{5};
ybar = varargin{6};
sx = varargin{7};
sy = varargin{8};
if nargin >= 9,
    cov_flag = 1;
    theta_cov = varargin{9};
    if nargin >=6,
        if varargin{10},
            cov_flag = 2; % hess flag
        end
    end
else
    cov_flag = 0;
    theta_cov = nan;
end
Ntheta = length(theta);
Nh = (Ntheta-Ny)/(Nx+Ny+1);
calllib('nnlib','nnlib_save_net',filename,Nx,Nh,Ny,theta,xbar,ybar,sx,sy,cov_flag,theta_cov');
out = {};


% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy,theta_cov] = nnlib('load_net',filename); 
%  (theta_cov returns empty if not in file)
% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy,theta_hess] = nnlib('load_net',filename,'theta_hess'); 
%  (theta_hess returns empty if not in file)
%  NOTE: theta_cov = inv(theta_hess)

function out = load_net(varargin)
filename = varargin{1};
if ~exist(filename,'file'),
    error('%s: load_net, file not found "%s"',mfilename,filename);
end

if length(varargin)>=2,
    theta_hess_flag = strcmpi('theta_hess',varargin{2});
else
    theta_hess_flag = false;
end

Nxptr = libpointer('uint32Ptr',0);
Nhptr = libpointer('uint32Ptr',0);
Nyptr = libpointer('uint32Ptr',0);
flagptr = libpointer('uint32Ptr',99999);
thetaptr = libpointer('doublePtr',nan);
xbarptr = libpointer('doublePtr',nan);
ybarptr = libpointer('doublePtr',nan);
sxptr = libpointer('doublePtr',nan);
syptr = libpointer('doublePtr',nan);
cov_thetaptr = libpointer('doublePtr',nan);
% call first time with flag=99999 to get sizes
calllib('nnlib','nnlib_load_net',filename,Nxptr,Nhptr,Nyptr,thetaptr,xbarptr,ybarptr,sxptr,syptr,flagptr,cov_thetaptr);
Nx = get(Nxptr,'value');
Nh = get(Nhptr,'value');
Ny = get(Nyptr,'value');
Ntheta = Ny+Ny*Nh+Nh+Nh*Nx;
flag = get(flagptr,'value');
thetaptr = libpointer('doublePtr',nan(Ntheta,1));
xbarptr = libpointer('doublePtr',nan(Nx,1));
ybarptr = libpointer('doublePtr',nan(Ny,1));
sxptr = libpointer('doublePtr',nan(Nx,1));
syptr = libpointer('doublePtr',nan(Ny,1));
if flag,
    cov_thetaptr = libpointer('doublePtr',nan(Ntheta,Ntheta));
else
    cov_thetaptr = libpointer('doublePtr',nan);
end
% call again to actually read data
calllib('nnlib','nnlib_load_net',filename,Nxptr,Nhptr,Nyptr,thetaptr,xbarptr,ybarptr,sxptr,syptr,flagptr,cov_thetaptr);

theta  = get(thetaptr,'value');
xbar = get(xbarptr,'value');
ybar = get(ybarptr,'value');
sx = get(sxptr,'value');
sy = get(syptr,'value');
if flag,
    cov_theta = get(cov_thetaptr,'value')';
    if xor(flag == 2,theta_hess_flag), % either cov_theta is cov_theta but want hess_theta, or cov_theta is hess_theta and want cov_theta
        cov_theta = inv(cov_theta);
    end
else
    cov_theta = [];
end
out = {Nx,Nh,Ny,theta,xbar,ybar,sx,sy,cov_theta};

function out = fit(varargin)

X = varargin{1};
Nh = varargin{2};
Y = varargin{3};
s = varargin{4};
MaxIter = varargin{5};
epsabs = varargin{6};
options = varargin(7:end);

[Nt,Nx] = size(X);
[Nt,Ny] = size(Y);

% determine flag based on size of s
if Ny==1,
    if numel(s)==1, % scalar, 1 x 1
        flag = 0;
    elseif max(size(s))==Nt, % Nt x 1 or Nt x 1 x 1
        flag = 2;
    else % 1 x Ny
        flag = 1;
    end
else
    flag = sum(size(s)>1);
end

if length(Nh)>1,
    theta = Nh;
    Ntheta = length(theta);
    Nh = calc_Nh(Nx,Ny,Ntheta); Nh = Nh{1};
    flag = flag+16; % start at initial theta
else
    Ntheta = calc_Ntheta(Nx,Nh,Ny); Ntheta = Ntheta{1};
    theta = nan(Ntheta,1);
end

for i = 1:length(options),
    switch(lower(options{i})),
        case 'do_cov',
            flag = flag-bitand(flag,132);
            flag = flag+4;
        case 'do_hess',
            flag = flag-bitand(flag,132);
            flag = flag+128;
        case 'bfgs',
            flag = flag+0;
        case 'fr',
            flag = flag+32;
        case 'pr',
            flag = flag+64;
        case 'nm',
            flag = flag+96;
        otherwise
            error('Unknown option "%s" to nnlib:fit',options{i});
    end
end

thetaPtr = libpointer('doublePtr',theta);
xbarptr = libpointer('doublePtr',nan(Nx,1));
ybarptr = libpointer('doublePtr',nan(Ny,1));
sxptr = libpointer('doublePtr',nan(Nx,1));
syptr = libpointer('doublePtr',nan(Ny,1));
cov_thetaPtr = libpointer('doublePtr',nan(Ntheta,Ntheta));

s = permute(s,ndims(s):-1:1);
calllib('nnlib','nnlib_fit',Nt,Nx,X',Nh,Ny,Y',s,flag,MaxIter,epsabs,thetaPtr,xbarptr,ybarptr,sxptr,syptr,cov_thetaPtr);

theta = get(thetaPtr,'value');
xbar = get(xbarptr,'value');
ybar = get(ybarptr,'value');
sx = get(sxptr,'value');
sy = get(syptr,'value');
theta_cov = get(cov_thetaPtr,'value')';
out = {theta,xbar,ybar,sx,sy,theta_cov};

function out = ell_debug(X,Nh,Y,s,theta,xbar,ybar,sx,sy)
[Nt,Nx] = size(X);
Ny = size(Y,2);
Ntheta = length(theta);

% determine flag based on size of s
if Ny==1,
    if numel(s)==1, % scalar, 1 x 1
        flag = 0;
    elseif max(size(s))==Nt, % Nt x 1 or Nt x 1 x 1
        flag = 2;
    else % 1 x Ny
        flag = 1;
    end
else
    flag = sum(size(s)>1);
end
s = permute(s,ndims(s):-1:1);

% double nnlib_ell_debug(const unsigned long int Nt, const unsigned long int Nx, const double *Z, 
% 			 const unsigned long int Nh, 
% 			 const unsigned long int Ny, const double *Y, 
% 			 const double *xbar, const double *ybar, 
% 			 const double *sx, const double *sy,
% 			 const double *s, const unsigned long int sflag,
% 			 const double *theta,
% 			 double *grad, double *hess);
gradPtr = libpointer('doublePtr',nan(Ntheta,1));
hessPtr = libpointer('doublePtr',nan(Ntheta,Ntheta));
ell = calllib('nnlib','nnlib_ell_debug',Nt,Nx,X',Nh,Ny,Y',xbar,ybar,sx,sy,s,flag,theta,gradPtr,hessPtr);
grad = get(gradPtr,'value');
hess = get(hessPtr,'value')';

out = {ell,grad,hess};
