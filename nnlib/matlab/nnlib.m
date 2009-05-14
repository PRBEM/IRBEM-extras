function varargout = nnlib(what,varargin)
%
% nnlib('save_training_set',filename,X,Y,s);
%
% [X,Y,s] = nnlib('load_training_set',filename);
%
% Y = nnlib('eval',X,Nh,theta,xbar,ybar,sx,sy)
% [Y,dY] = nnlib('eval',X,Nh,theta,xbar,ybar,sx,sy,dY_flag,theta_cov)
% Y = nnlib('eval',X,net))
% [Y,dY] = nnlib('eval',X,net,dY_flag)
%    dY_flag & 3 = 0: don't compute dY (theta_cov ignored)
%    dY_flag & 3 = 1: for 1 x Ny dY (standard errors)
%    dY_flag & 3 = 2: for Ny x Ny dY (covariance matrix)
%    dY_flag & 4 = 4: theta_cov is theta_hess (not need in net syntax)
% [ell,grad,hess] = nnlib('ell_debug',X,Nh,Y,s,theta,xbar,ybar,sx,sy)
%
% nnlib('save_net',filename,Nx,Ny,theta,xbar,ybar,sx,sy,theta_cov);
% nnlib('save_net',filename,Nx,Ny,theta,xbar,ybar,sx,sy,theta_hess,'hess_theta');
%
% net = nnlib('load_net',filename); 
% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy] = nnlib('load_net',filename); 
% (theta_cov not read, even if present in file)
% net = nnlib('load_net',filename,'theta_cov'); 
% [Nx,Nh,Ny,theta,xbar,ybar,sx,sy,theta_cov] = nnlib('load_net',filename); 
%  (theta_cov returns empty if not in file)
% net = nnlib('load_net',filename,'theta_hess'); 
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
    if ispc,
        libfile = which('nnlib.dll');
    else
        libfile = which('nnlib.so');
    end
    hfile = which('nnlib.h');
    loadlibrary(libfile,hfile,'alias','nnlib');
end

varargout = cell(1,nargout);

switch(lower(what)),
    case {'save_training_set'},
        [varargout{:}] = save_training_set(varargin{:});
    case {'load_training_set'},
        [varargout{:}] = load_training_set(varargin{:});
    case {'fit'},
        [varargout{:}] = fit(varargin{:});        
    case {'ell_debug'},
        [varargout{:}] = ell_debug(varargin{:});        
    case {'save_net'},
        [varargout{:}] = save_net(varargin{:});
    case {'load_net'},
        [varargout{:}] = load_net(varargin{:});
    case {'eval'},
        [varargout{:}] = eval(varargin{:});
    case {'ntheta'},
        [varargout{:}] = calc_Ntheta(varargin{:});
    case {'nh'},
        [varargout{:}] = calc_Nh(varargin{:});
    otherwise
        error('%s: Unknown "what": "%s"',mfilename,what);
end

function Ntheta = calc_Ntheta(Nx,Nh,Ny)
Ntheta = Ny+Ny*Nh+Nh+Nh*Nx;

function Nh = calc_Nh(Nx,Ny,Ntheta)
Nh = (Ntheta-Ny)/(Ny+Nx+1);

function save_training_set(filename,X,Y,s)

Nt = size(X,1);
Nx = size(X,2);
Ny = size(Y,2);

if Ny==1,
    if numel(s)==1, % scalar, 1 x 1
        flag = 0;
    elseif numel(s) == Nt, % Nt x 1 or Nt x 1 x 1
        flag = 2;
    elseif numel(s) == Ny,
        flag = 1;
    else
        error('Unable to deduce meaning of s with %d elements',numel(s));
    end
else
    flag = sum(size(s)>1);
end
s = permute(s,ndims(s):-1:1);
calllib('nnlib','nnlib_save_training_set',filename,Nt,Nx,X',Ny,Y',s,flag);

function [X,Y,s] = load_training_set(varargin)
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
Nx = Nxptr.value;
Ny = Nyptr.value;
Nt = Ntptr.value;
flag = flagptr.value;
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

X = Xptr.value';
Y = Yptr.value';
s = reshape(sptr.value,size_s); % libpointers can only hold matrices, no 3-D stuff
s = permute(s,ndims(s):-1:1);

function [Y,dY] = eval(X,varargin)
Nt = size(X,1);
Nx = size(X,2);
if length(varargin) > 2,
    Nh = varargin{1};
    theta = varargin{2};
    xbar = varargin{3};
    ybar = varargin{4};
    sx = varargin{5};
    sy = varargin{6};
    dY_flag = varargin{7};
    theta_cov = varargin{8};
else
    if length(varargin)==2,
        dY_flag = varargin{2};
    else
        dY_flag = 0;
    end
    net = varargin{1};
    Nh = net.Nh;
    theta = net.theta;
    xbar = net.xbar;
    ybar = net.ybar;
    sx = net.sx;
    sy = net.sy;
    if isfield(net,'theta_cov'),
        theta_cov = net.theta_cov;
    end
    if isfield(net,'theta_hess'),
        theta_hess = net.theta_hess;
        theta_cov = theta_hess;
        dY_flag = bitor(dY_flag,4); % turn on flag indicating theta_cov is theta_hess
    end
end

Ntheta = length(theta);
Ny = (Ntheta-(Nx+1)*Nh)/(Nh+1);

Yptr = libpointer('doublePtr',nan(Ny,Nt));
if nargin>=8,
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

Y = Yptr.value';
if nargout>=2,
    if bitand(dY_flag,3),
        dY = dYptr.value';
        if bitand(dY_flag,3) == 2,
            dY = reshape(dY,[Nt,Ny,Ny]);
        end
    else
        dY = [];
    end
end

function save_net(filename,Nx,Ny,theta,xbar,ybar,sx,sy,theta_cov,hess_theta_flag)
if nargin >= 9,
    if nargin >=10,
        cov_flag = 1+strcmpi(hess_theta_flag,'hess_theta'); % hess flag
    else
        cov_flag = 1;
    end
else
    cov_flag = 0;
    theta_cov = nan;
end
Ntheta = length(theta);
Nh = (Ntheta-Ny)/(Nx+Ny+1);
calllib('nnlib','nnlib_save_net',filename,Nx,Nh,Ny,theta,xbar,ybar,sx,sy,cov_flag,theta_cov');

function varargout = load_net(filename,theta_hess_flag)
do_struct = (nargout == 1);
cov_name = 'theta_cov';
if nargin>=2,
    do_cov = true;
    theta_hess_flag = strcmpi('theta_hess',theta_hess_flag);
else
    do_cov = (nargout > 8);
    theta_hess_flag = false;
end

if theta_hess_flag,
    cov_name = 'theta_hess';
end

if ~exist(filename,'file'),
    error('%s: load_net, file not found "%s"',mfilename,filename);
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
Nx = Nxptr.value;
Nh = Nhptr.value;
Ny = Nyptr.value;
Ntheta = Ny+Ny*Nh+Nh+Nh*Nx;
if do_cov,
    flag = flagptr.value;
else % don't load flag
    flag = 0;
    flagptr.value = flag;
end
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

theta  = thetaptr.value;
xbar = xbarptr.value;
ybar = ybarptr.value;
sx = sxptr.value;
sy = syptr.value;
if flag,
    cov_theta = cov_thetaptr.value';
    if xor(flag == 2,theta_hess_flag), % either cov_theta is cov_theta but want hess_theta, or cov_theta is hess_theta and want cov_theta
        cov_theta = inv(cov_theta);
    end
else
    cov_theta = [];
end

if do_struct,
    varargout = {struct('Nx',Nx,'Nh',Nh,'Ny',Ny,'theta',theta,'xbar',xbar,'ybar',ybar,'sx',sx,'sy',sy)};
    if do_cov,
        varargout{1}.(cov_name) = cov_theta;
    end
else
    varargout = {Nx,Nh,Ny,theta,xbar,ybar,sx,sy};
    if do_cov,
        varargout{end+1} = cov_theta;
    end
end


function [theta,xbar,ybar,sx,sy,theta_cov] = fit(X,Nh,Y,s,MaxIter,epsabs,varargin)

options = varargin;

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
    Nh = calc_Nh(Nx,Ny,Ntheta);
    flag = flag+16; % start at initial theta
else
    Ntheta = calc_Ntheta(Nx,Nh,Ny);
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

theta = thetaPtr.value;
xbar = xbarptr.value;
ybar = ybarptr.value;
sx = sxptr.value;
sy = syptr.value;
theta_cov = cov_thetaPtr.value';

function [ell,grad,hess] = ell_debug(X,Nh,Y,s,theta,xbar,ybar,sx,sy)
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
if nargout >=2,
    gradPtr = libpointer('doublePtr',nan(Ntheta,1));
else
    gradPtr = libpointer('doublePtr'); % null pointer, ignore
end
if nargout >= 3,
    hessPtr = libpointer('doublePtr',nan(Ntheta,Ntheta));
else
    hessPtr = libpointer('doublePtr'); % null pointer, ignore
end
ell = calllib('nnlib','nnlib_ell_debug',Nt,Nx,X',Nh,Ny,Y',xbar,ybar,sx,sy,s,flag,theta,gradPtr,hessPtr);
grad = get(gradPtr,'value');
hess = get(hessPtr,'value')';
