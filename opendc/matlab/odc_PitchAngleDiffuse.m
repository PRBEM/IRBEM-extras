function [psd,times,alpha0] = odc_PitchAngleDiffuse(alpha0,f0,Daa,times)
% Demonstrate decay of a pitch-angle distribution
% psd = odc_PitchAngleDiffuse(alpha0,f0,Daa,times)
% pitch-angle diffusion eigenmodes
% Assuming elastic pitch-angle-diffusion
% From Schulz and Lanzerotti eq 2.16, page 56
% df/dt = d[x T(y) Dxx df/dx]/dx/T(y)/x
% x = cos(alpha0), y = sin(alpha0) = sqrt(1-x^2)
% alpha0 - equatorial pitch angle grid, degrees
% (starting at alpha0=0 is bad as it introduces strong dependence on grid
%  spacenig. Use a finite loss cone angle instead)
% f0 - function f0(alpha0) is initial condition (in arbitrary units)
% Daa - function Daa(alpha0) is diffusion coefficient (1/time))
% times- array of times to spit out (same time unit as Daa)
% psd is length(times) x length(alpha0)
% [psd,times,alpha0] = ... also returns grids
% use [] to supply defaults

if (nargin < 1) || isempty(alpha0),
    alpha0 = 5:90;
end
alpha0 = alpha0(:)'; % row vector

if (nargin < 2) || isempty(f0),
    f0 = @(X)ones(size(X));
end

if (nargin < 3) || isempty(Daa),
    Daa = @(a)ones(size(a)); % constant
end

if (nargin < 4) || isempty(times),
    times = 0:0.01:1;
end

Dxx = @(x)Daa(acosd(x)).*(1-x.^2); % Dxx = Daa*y^2

%C(x,t,u,Du/Dx) * Du/Dt = x^(-M) * D(x^M * F(x,t,u,Du/Dx))/Dx + S(x,t,u,Du/Dx)

% M = 1
% C = T(y)
% F = T(y) Dxx DuDx
% S = 0

[XMESH,isort] = sort(cosd(alpha0));

M = 1;
psd = pdepe(M,@(X,T,U,DUDX)PDEFUN(X,T,U,DUDX,Dxx),...
    @(x)f0(acosd(x)),@(XL,UL,XR,UR,T)BCFUN(XL,UL,XR,UR,T,Dxx),...
    XMESH,times);
psd(:,isort) = psd; % undo sort by x

function [PL,QL,PR,QR] = BCFUN(XL,UL,XR,UR,T,Dxx)
% P(x,t,u) + Q(x,t) * F(x,t,u,Du/Dx) = 0

% du/dx(x=0) = 0
PL = 0;
QL = (Dxx(XL).*SL_T(sqrt(1-XL.^2))).^-1; % cancel out T Dxx in F

% u(x=1) = 0; loss cone
PR = UR;
QR = 0;


function [C,F,S] = PDEFUN(X,T,U,DUDX,Dxx)
T = SL_T(sqrt(1-X.^2)); % T(y);
C = diag(T);
F = T.*Dxx(X).*DUDX;
S = zeros(size(X));

function T = SL_T(y)
% Schulz & Lanzerotti's T(y)
% Schulz & Lanzerotti constants
T0 = 1+1/(2*sqrt(3))*log(2+sqrt(3)); % S&L 1.28a
T1 = pi/6*sqrt(2); % S&L 1.28b
T = T0-0.5*(T0-T1)*(y+sqrt(y)); % 1/4 bounce integral of 1
