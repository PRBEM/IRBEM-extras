function inst_info = rfl_init_Eseparable(inst_info)
% inst_info = rfl_init_Eseparable(inst_info)
% populate methods and extra properties for
% an [E],[...] response type
% inst_info is assumed to be a fully-initialized inseparable R(E,theta,phi) sensor
% this function just overloads the definition of R

% separable variables overload all methods because
% R is different, and the integral weights methods
% can be done much faster when integrals are separated

inst_info.internal.R = @(inst_info,E,theta,phi)inst_info.internal.RE(inst_info,E).*inst_info.internal.A(inst_info,theta,phi);

inst_info.make_hEthetaphi = @make_hEthetaphi_Eseparable;
inst_info.make_hEtheta = @make_hEtheta_Eseparable;
inst_info.make_hE = @make_hE_Eseparable;
inst_info.make_hthetaphi = @make_hthetaphi_Eseparable;
inst_info.make_htheta = @make_htheta_Eseparable;
inst_info.make_hEalphabeta = @make_hEalphabeta_Eseparable;
inst_info.make_hEalpha = @make_hEalpha_Eseparable;
inst_info.make_hEiso = @make_hEiso_Eseparable;
inst_info.make_halphabeta = @make_halphabeta_Eseparable;
inst_info.make_halpha = @make_halpha_Eseparable;

% separables have a new structure "internal"
% with (protected) methods:
% inst_info.internal.hAthetaphi(inst_info,thetagrid,phigrid,options), used for hE*hA
inst_info.internal.hAthetaphi = @make_hAthetaphi;
% inst_info.internal.hAtheta(inst_info,thetagrid,options), used for hE*hA
inst_info.internal.hAtheta = @make_hAtheta;
% inst_info.internal.hAalphabeta(inst_info,...,options), used for hE*hA
inst_info.internal.hAalphabeta = @make_hAalphabeta;
% inst_info.internal.hAalpha(inst_info,...,options), used for hE*hA
inst_info.internal.hAalpha = @make_hAalpha;
% inst_info.internal.hA0(inst_info,options) for isotropic flux (G)
inst_info.internal.hA0 = @make_hA0; % - will become a constant of type double before return 
inst_info.internal.merge_hE_hangles = @merge_hE_hangles; % useful for derived classes to call, too
% sensor-specific protected functions/methods:
% inst_info.internal.A(inst_info,theta,phi) - area function, used for R
% inst_info.internal.RE(inst_info,E) - eps(E), used for R
% inst_info.internal.hE0 - hE for flat spectrum (1 for integral channels)
% inst_info.internal.hE(inst_info,Egrid,options), used for hE*hA


% these init functions define internal.RE, internal.hE, and internal.hE0
switch(inst_info.E_TYPE),
    case 'TBL',
        inst_info = rfl_init_TBL(inst_info);
    case 'INT',
        inst_info = rfl_init_INT(inst_info);
    case 'WIDE',
        inst_info = rfl_init_WIDE(inst_info);
    case 'DIFF',
        inst_info = rfl_init_DIFF(inst_info);
    otherwise
        error('E_TYPE %s not defined yet',inst_info.E_TYPE);
end

% these functions define inst_info.internal.A, .hA* and .hE*
% and may overload inst_info.make_hEalphabeta, etc.

switch(inst_info.RESP_TYPE),
    case '[E],[TH,PH]',
        switch(inst_info.TP_TYPE),
            case 'TBL',
                inst_info.internal.A = @(inst_info,theta,phi)interpn(inst_info.TH_GRID,inst_info.PH_GRID,inst_info.A,theta,phi,'linear',0); % zero outside grid                
            case 'RECT_TELE',
                inst_info = init_telescope_rect(inst_info);
            otherwise
                error('TP_TYPE %s not defined yet',inst_info.TP_TYPE);
        end
    case {'[E],[TH]','[E]'},
        inst_info = rfl_init_Eseparable_csym(inst_info); % overload for cylindrically symmetric case
    otherwise
        error('RESP_TYPE %s not defined yet',inst_info.RESP_TYPE);
end

if isa(inst_info.internal.hA0,'function_handle'), % this is a constant, once initialized
    inst_info.internal.hA0 = inst_info.internal.hA0(inst_info);
end

function inst_info = init_telescope_rect(inst_info)
% initialize rectangular 2-element telescope
% overloads inst_info.internal.A, sets some inst_info.internal fields, too
% W1,H1,W2,H2, D, BIDIRECTIONAL

% Using Sullivan's 1971 paper as updated in ~2010
inst_info.internal.X = @(zeta,a1,a2) min(zeta+a1/2,a2/2) - max(zeta-a1/2,-a2/2);
inst_info.internal.A = @telescope_rect_A;
D = inst_info.D; % Sullivan's ell
alpha = (inst_info.H1 + inst_info.H2)/2;
beta = (inst_info.W1 + inst_info.W2)/2;
gamma = (inst_info.H1 - inst_info.H2)/2;
delta = (inst_info.W1 - inst_info.W2)/2;

% Sullivan's (11)
inst_info.internal.hA0 = D^2*log(((D^2+alpha^2+delta^2)/(D^2+alpha^2+beta^2))*(D^2+gamma^2+beta^2)/(D^2+gamma^2+delta^2)) + ...
     2*alpha*sqrt(D^2+beta^2) *atan(alpha/sqrt(D^2+beta^2)) +2*beta *sqrt(D^2+alpha^2)*atan(beta /sqrt(D^2+alpha^2)) + ...
    -2*alpha*sqrt(D^2+delta^2)*atan(alpha/sqrt(D^2+delta^2))-2*beta *sqrt(D^2+gamma^2)*atan(beta /sqrt(D^2+gamma^2)) + ...
    -2*gamma*sqrt(D^2+beta^2) *atan(gamma/sqrt(D^2+beta^2)) -2*delta*sqrt(D^2+alpha^2)*atan(delta/sqrt(D^2+alpha^2)) + ...
     2*gamma*sqrt(D^2+delta^2)*atan(gamma/sqrt(D^2+delta^2))+2*delta*sqrt(D^2+gamma^2)*atan(delta/sqrt(D^2+gamma^2));
if inst_info.internal.bidirectional,
    inst_info.internal.hA0 = inst_info.internal.hA0*2;
end

function A = telescope_rect_A(inst_info,theta,phi)
% effective area for 2-element rectangular telescope
% right out of Sullivan's updated paper, equation (13)
A = zeros(size(theta));

tantheta = tand(theta);
zeta = inst_info.D*tantheta.*cosd(phi);
eta = inst_info.D*tantheta.*sind(phi);

X = inst_info.internal.X(zeta,inst_info.W1,inst_info.W2);
Y = inst_info.internal.X(eta,inst_info.H1,inst_info.H2);
if inst_info.internal.bidirectional,
    theta = min(theta,180-theta);
    costheta = cosd(theta);
else
    costheta = cosd(min(90,theta));
end
f = find((costheta>0) & (X>0) & (Y>0)); % apply Heaviside implicitly and resolve tan(90)=inf
if any(f(:)),
    A(f) = costheta(f).*X(f).*Y(f); % eq 13 w/o Heaviside functions
end

function [hAthetaphi,result_code] = make_hAthetaphi(inst_info,thetagrid,phigrid,options)
dcost = rfl_make_deltas(-cosd(thetagrid),options);
dp = rfl_make_deltas(phigrid,options)*pi/180;
[dcost,dp] = ndgrid(dcost,dp);
[theta,phi] = ndgrid(thetagrid,phigrid);

result_code = 1; % success
A = inst_info.internal.A(inst_info,theta,phi);
h = A.*dcost.*dp;
hAthetaphi = h; % success, returns h

function [hAtheta,result_code] = make_hAtheta(inst_info,thetagrid,options)
phigrid = rfl_make_grid(0,360,'phi',options);
[hAthetaphi,result_code] = inst_info.internal.hAthetaphi(inst_info,thetagrid,phigrid,options);
if result_code == 1,
    hAtheta = sum(hAthetaphi,2);
else
    hAtheta = nan;
end

function [hA0,result_code] = make_hA0(inst_info)
thetagrid = rfl_make_grid(0,180,'theta',[]);
[hAtheta,result_code] = inst_info.internal.hAtheta(inst_info,thetagrid,[]);
if result_code == 1,
    hA0 = sum(hAtheta);
else
    hA0 = nan;
end

function [hAalphabeta,result_code] = make_hAalphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
if isfield(options,'acute'),
    acute = options.acute;
else
    acute = alphagrid(end)<=90;
end
dcosa = rfl_make_deltas(-cosd(alphagrid),options);
db = rfl_make_deltas(betagrid,options)*pi/180;
dt = rfl_make_deltas(tgrid,options);
[dcosa,db] = ndgrid(dcosa,db);
[alpha,beta] = ndgrid(alphagrid,betagrid);
theta = nan(size(alpha));
phi = theta;
% prepare to loop
result_code = 1; % success
hAalphabeta = nan; % on error, returns a nan
h = 0;
for it = 1:length(tgrid),
    [theta(:),phi(:),result_code] = rfl_alphabeta2thetaphi(alpha(:),beta(:),alpha0(it),beta0(it),phib(it));
    if result_code ~= 1,
        return
    end
    h = h+inst_info.internal.A(inst_info,theta,phi)*dt(it);
    if acute, % add pi-alpha point as well
        [theta(:),phi(:),result_code] = rfl_alphabeta2thetaphi(180-alpha(:),beta(:),alpha0(it),beta0(it),phib(it));
        if result_code ~= 1,
            return
        end
        h = h+inst_info.internal.A(inst_info,theta,phi)*dt(it);
    end
end
h = h.*dcosa.*db;
hAalphabeta = h; % success, returns h

function [hAalpha,result_code] = make_hAalpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options)
betagrid = rfl_make_grid(0,360,'beta',options);
[hAalphabeta,result_code] = inst_info.internal.hAalphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);
if result_code == 1,
    hAalpha = sum(hAalphabeta,2);
else
    hAalpha = nan;
end

function [hEthetaphi,result_code] = make_hEthetaphi_Eseparable(inst_info,Egrid,thetagrid,phigrid,options)
% combines results of make_hE and make_thetaphi
result_code = 1;
hE = inst_info.internal.hE(inst_info,Egrid,options);
hAthetaphi = inst_info.internal.hAthetaphi(inst_info,thetagrid,phigrid,options);

hEthetaphi = merge_hE_hangles(hE,hAthetaphi);

function [hthetaphi,result_code] = make_hthetaphi_Eseparable(inst_info,thetagrid,phigrid,options)
result_code = 1;
hAthetaphi = inst_info.internal.hAthetaphi(inst_info,thetagrid,phigrid,options);
hthetaphi = hAthetaphi*inst_info.internal.hE0;% includes XCAL = CROSSCALIB

function [hEtheta,result_code] = make_hEtheta_Eseparable(inst_info,Egrid,thetagrid,options)
% combines results of make_hE and make_theta
result_code = 1;
hE = inst_info.internal.hE(inst_info,Egrid,options);
hAtheta = inst_info.internal.hAtheta(inst_info,thetagrid,options);
hEtheta = merge_hE_hangles(hE,hAtheta);

function [htheta,result_code] = make_htheta_Eseparable(inst_info,thetagrid,options)
% combines results of make_hE and make_theta
result_code = 1;
hAtheta = inst_info.internal.hAtheta(inst_info,thetagrid,options);
htheta = hAtheta*inst_info.internal.hE0;% includes XCAL = CROSSCALIB

function [hE,result_code] = make_hE_Eseparable(inst_info,Egrid,options)
[hE,result_code] = inst_info.internal.hE(inst_info,Egrid,options);
hE = hE*inst_info.internal.hA0;

function [hEalphabeta,result_code] = make_hEalphabeta_Eseparable(inst_info,Egrid,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
% combines results of make_hE and make_alphabeta
result_code = 1;
hE = inst_info.internal.hE(inst_info,Egrid,options);
hAalphabeta = inst_info.internal.hAalphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);
hEalphabeta = merge_hE_hangles(hE,hAalphabeta);

function [halphabeta,result_code] = make_halphabeta_Eseparable(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options)
result_code = 1;
hAalphabeta = inst_info.internal.hAalphabeta(inst_info,alphagrid,betagrid,tgrid,alpha0,beta0,phib,options);
halphabeta = hAalphabeta*inst_info.internal.hE0; % includes XCAL = CROSSCALIB

function [hEalpha,result_code] = make_hEalpha_Eseparable(inst_info,Egrid,alphagrid,tgrid,alpha0,beta0,phib,options)
% combines results of make_hE and make_alpha
result_code = 1;
hE = inst_info.internal.hE(inst_info,Egrid,options);
hAalpha = inst_info.internal.hAalpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options);
hEalpha = merge_hE_hangles(hE,hAalpha);

function [halpha,result_code] = make_halpha_Eseparable(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options)
result_code = 1;
hAalpha = inst_info.internal.hAalpha(inst_info,alphagrid,tgrid,alpha0,beta0,phib,options);
halpha = hAalpha*inst_info.internal.hE0; % includes XCAL = CROSSCALIB

function [hEiso,result_code] = make_hEiso_Eseparable(inst_info,Egrid,tgrid,options)
% combines results of make_hE and hA0
[hE,result_code] = inst_info.internal.hE(inst_info,Egrid,options);
if result_code == 1,
    dt = rfl_make_deltas(tgrid,options);
    hEiso = hE*sum(dt);
else
    hEiso = nan;
end

function hEthetaphi = merge_hE_hangles(hE,hangles)

% hE has shape: NE x 1
% hangles has shape: N1 x N2
% hEthetaphi needs shape NE x N1 x N2

NE = length(hE);
[N1,N2] = size(hangles);

hangles = shiftdim(hangles,-1); % introduce singleton first dimension
hEthetaphi = repmat(hE,[1 N1 N2]).*repmat(hangles,[NE,1,1]);

%%%%%%%%%%%%%%% energy response functions %%%%%%%%%%%%%%%%%%%%%%%%

function inst_info = rfl_init_DIFF(inst_info)
% idealized differential energy channel

inst_info.internal.RE = @(inst_info,E)(E==inst_info.E0)*inst_info.EPS; % DE gets ignored, which is bad, but unaviodable
inst_info.internal.hE = @rfl_hE_diff;
inst_info.internal.hE0 = inst_info.DE*inst_info.EPS/inst_info.CROSSCALIB; % hE for flat spectrum

function [hE,result_code] = rfl_hE_diff(inst_info,Egrid,options)
% options is ignored
result_code = 1;
NE = length(Egrid);
hE = zeros(size(Egrid));
I = rfl_get_list_neighbors(Egrid,inst_info.E0);

% Egrid(I(1)) <= inst_info.E0 < Egrid(I(2))
% or I(2) = 0 or both are 0

i = I(1);
if (i>=1) && (i < NE), % E(i) <= E0 < E(i+1)
    hE(i) = inst_info.DE*inst_info.EPS*(Egrid(i+1) - inst_info.E0) / (Egrid(i+1)-Egrid(i)) ;
end

i = I(2);
if (i>=2) && (i <= NE), % E(i-1) < E0 < E(i)
    hE(i) = inst_info.DE*inst_info.EPS*(inst_info.E0-Egrid(i-1)) / (Egrid(i)-Egrid(i-1));
end

hE = hE / inst_info.CROSSCALIB;

function inst_info = rfl_init_INT(inst_info)
% idealized integral energy channel

inst_info.internal.RE = @(inst_info,E)(E>=inst_info.E0)*inst_info.EPS;
inst_info.internal.hE = @rfl_hE_int;
inst_info.internal.hE0 = inst_info.EPS/inst_info.CROSSCALIB; % should be inf for flat spectrum, but assume energy bandwidth==1

function [hE,result_code] = rfl_hE_int(inst_info,Egrid,options)
result_code = 1;
NE = length(Egrid);
I = rfl_get_list_neighbors(Egrid,inst_info.E0);

% Egrid(I(1)) <= inst_info.E0 < Egrid(I(2))
% or I(2) = 0 or both I=0

dE = rfl_make_deltas(Egrid,options);
hE = dE; % default, eventually only kept for Egrid>Egrid(I(2))
hE(1:(I(1))) = 0;

i = I(1);
if (i>=1) && (i < NE), % E(i) <= E0 < E(i+1)
    hE(i) = (Egrid(i+1)-inst_info.E0)^2/(2*(Egrid(i+1)-Egrid(i)));
end

i = I(2);
if (i>=2) && (i <= NE), % E(i-1) < E0 < E(i)
    hE(i) = dE(i) - (inst_info.E0-Egrid(i-1))^2/(2*(Egrid(i)-Egrid(i-1)));
end
hE = hE * inst_info.EPS / inst_info.CROSSCALIB;


function inst_info = rfl_init_WIDE(inst_info)
% idealized wide differential energy channel

inst_info.internal.RE = @(inst_info,E)((E>=inst_info.E0)&(E<=inst_info.E1))*inst_info.EPS;
inst_info.internal.hE = @rfl_hE_wide;
inst_info.internal.hE0 = (inst_info.E1-inst_info.E0)*inst_info.EPS/inst_info.CROSSCALIB; % hE for flat spectrum

function [hE,result_code] = rfl_hE_wide(inst_info,Egrid,options)
% treat as difference between two integral channels
% at E1 and E0
result_code = 1;
resptmp = inst_info;
resptmp.E0 = resptmp.E1;
hE = rfl_hE_int(inst_info,Egrid,options)-rfl_hE_int(resptmp,Egrid,options);

function inst_info = rfl_init_TBL(inst_info)
% idealized wide differential energy channel

% min/max forces plateau extrapolation in E
inst_info.internal.RE = @(inst_info,E)interp1(inst_info.E_GRID,inst_info.EPS,min(max(E,inst_info.E_GRID(1)),inst_info.E_GRID(end)),'linear',0); % set to zero outside grid

inst_info.internal.hE = @rfl_hE_tbl;
inst_info.internal.hE0 = sum(inst_info.E_GRID(:).*inst_info.EPS(:).*rfl_make_deltas(inst_info.E_GRID(:)))/inst_info.CROSSCALIB; % hE for flat spectrum

function [hE,result_code] = rfl_hE_tbl(inst_info,Egrid,options)
result_code = 1;
dE = rfl_make_deltas(Egrid,options);
hE = inst_info.internal.RE(inst_info,Egrid).*dE/inst_info.CROSSCALIB;

