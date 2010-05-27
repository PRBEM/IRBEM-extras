function varargout = invlib(what,varargin)
% invlib provides access to the irbem/invlib library
% the following routines are supported: ana_spec_inv, omni2uni
%
% ana_spec_inv:
% fit = invlib('ana_spec_inv',y,dy,Egrid,G,dt,b,Eout,...)
% calls ana_spec_inv_multi as defined in invlib.pdf
% y and b are NT x NY, dt is NT x 1
% dy is NY x 1 
% Egrid is NE x 1, Eout is NEout x 1
% G is defined as in int dE j(E) G(E) dE, size NY x NE, 
% options:
% select minimizer : 'BFGS' (default), 'FR', 'PR', 'NM'
% select analytical spectrum functions: 'EXP','PL','PLE','RM','RM2' 
%   (at least one function required)
% set 'trapz' to use trapezoidal energy integral (default)
% set 'plateau' to use plateau weights in energy integral
% set 'G=GdE' to assume dE included in G
% set rest_energy for RM and RM2: ...,'rest_energy',rest_energy,...
% set Ebreak for PLE: ...,'Ebreak',Ebreak,...
% set E0 for PLE: ...,'E0',E0,...
% set maximum fit iterations: ...,'MaxIter',MaxIter,...
% set verbose output text file: ...,'outfile',filename,...
% set 'append' to output file: ...,'append',...
%
% fit has fields:
% flux, dlogflux, lambda, result, result_codes as returned by ana_spec_inv_multi
% flux.exp, flux.pl, flux.ple, flux.rm, flux.rm2 have the following fields:
% flux.exp.q(t,m) m'th fit parameter for t'th sample
% flux.exp.Nq = # of free parameters in analytical fit
% flux.exp.flux(q(t,:),E) returns fit flux at energies E for sample t
% flux.exp.ell = negative log likelihood of fit
% flux.exp.weight = weight of fit in combined solution (i.e., in fit.flux)
% flux.exp.lambda(t,i) = expected counts in t'th sample, i'th channel for fit
%
% pc_spec_inv: perform principal component spectral inverse
% fit = invlib('pc_spec_inv',y,dy,Egrid,G,dt,b,mean_log_flux,basis_vectors,basis_variance...)
% common arguments have same meaning as for ana_spec_inv
% mean_log_flux - NE x 1 vector of mean log flux
% basis_vectors - NE x Nq matrix of basis vectors (columns)
% basis_variance - Nq x 1 vector of variance in log flux explained by each
%    basis vector
% options: supports appropriate options for ana_spec_inv, plus
% 'num_bases',num_bases - select a reduced number of bases to use
% fit has fields:
% flux, dlogflux, lambda, result, result_codes as returned by pc_spec_inv_multi
% fit.ell = negative log likelihood of fit
% flux.q(t,m) m'th PC coefficient for t'th sample
%
% omni2uni:
% [uniflux,dloguniflux,result] = invlib('omni2uni',omniflux,dlogomniflux,method,...)
% calls omni2uni as defined in invlib.pdf
% (note, if omniflux is a vector, then it'll invert each omniflux)
% omniflux - Estimated isotropic unidirectional flux 
%  (for TEM-1 method, use #/cm^2/sr/s/keV)
% dlogomniflux - relative error for omniflux (dimensionless)
% method is one of: (required options given in parentheses)
%  'TEM1' - uses TEM1 analytical climatology for electrons (L, E, B/B0)
%  'Vampola' - uses Vampola's L-dependent sin^n for electrons (L)
% key,value options (default given in parens):
%  'keV' - Energy of particle flux, keV
%  'MeV' - Energy of particle flux, MeV
%  'Lm' - McIlwain L of locally mirroring particles, in OPQ
%  'B/B0' - B/B0 of spacecraft location, in OPQ
%  'NA' - number of grid points for angular integral (50)
%  'MaxIter' - maximum number of iterations for minimizer (1,000)
% keyword options for choice of root finder
%  'PHS' - Powell's Hybrid method, scaled (default)
%  'PH' - Powell's Hybrid method
%  'PHNHS' - Powell's Hybrid method w/ numerical Hessian, scaled 
%  'PHNH' - Powell's Hybrid method w/ numerical Hessian
%  'Newton' - Newton's method
%  'GN' - Pseudo-global Newton's method
%  'DN' - Discrete Newton's method
%  'Broyden' - Broyden's algorithm, numerical Hessian
% uniflux - locally-mirroring unidirectional flux (e.g., in #/keV/cm^2/sr/s)
% dloguniflux - standard error of natural log of uniflux (dimensionless)
% result - return codem fro library call
%  1 - success, 0 or negative is an error code (see invlib.pdf)

load_invlib;

varargout = cell(1,nargout);

switch(lower(what)),
    case {'ana_spec_inv'},
        [varargout{:}] = ana_spec_inv(varargin{:});
    case {'pc_spec_inv'},
        [varargout{:}] = pc_spec_inv(varargin{:});
    case {'omni2uni'},
        [varargout{:}] = omni2uni(varargin{:});
    otherwise
        error('%s: Unknown "what": "%s"',mfilename,what);
end

function [NT,NY,NE] = spec_inv_check_args(y,dy,Egrid,H0,dt,b)
[NT,NY] = size(y);
NE = length(Egrid);
if numel(dy) ~= NY,
    error('dy must have size NYx1');
end
if (numel(H0) ~= NY*NE) || (size(H0,1) ~= NY),
    error('H0 must have size NY x NE');
end
if numel(dt) ~= NT,
    error('dt must have size NT x 1');
end
if (numel(b) ~= NT*NY) || (size(b,1) ~= NT),
    error('b must have size NT x NY');
end

function options = parse_options(varargin)
% constants from specinv.h
options.ASI_FXN_PL   = (1);
options.ASI_FXN_EXP  = (2);
options.ASI_FXN_RM   = (4);
options.ASI_FXN_PLE  = (8);
options.ASI_FXN_RM2  = (16);
options.ASI_MAX_POW2 = 4;
options.ASI_MAX_NQ = 10;

% constants from optim.h
options.OPTIM_MIN_BFGS  = (0);
options.OPTIM_MIN_FR    = (1);
options.OPTIM_MIN_PR    = (2);
options.OPTIM_MIN_NM    = (3);

options.num_bases = inf;
options.minimizer = 0; % BFGS
options.MaxIter = 1000;
options.fxn_bitmap = 0;
options.outfile = '';
options.rest_energy = 0.511; % electron rest energy, MeV
options.Ebreak = 100; % power-law to exponential tail transition energy, MeV
options.E0 = 345; % power-law to exponential tail transition energy, MeV
options.append = false; % append to text file?
options.dE_mode = 1; % default - trapz


i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'bfgs',
            options.minimizer = options.OPTIM_MIN_BFGS;
        case 'fr',
            options.minimizer = options.OPTIM_MIN_FR;
        case 'pr',
            options.minimizer = options.OPTIM_MIN_PR;
        case 'nm',
            options.minimizer = options.OPTIM_MIN_NM;
        case 'pl',
            options.fxn_bitmap = bitor(options.fxn_bitmap,options.ASI_FXN_PL);
        case 'exp',
            options.fxn_bitmap = bitor(options.fxn_bitmap,options.ASI_FXN_EXP);
        case 'rm',
            options.fxn_bitmap = bitor(options.fxn_bitmap,options.ASI_FXN_RM);
        case 'ple',
            options.fxn_bitmap = bitor(options.fxn_bitmap,options.ASI_FXN_PLE);
        case 'rm2',
            options.fxn_bitmap = bitor(options.fxn_bitmap,options.ASI_FXN_RM2);
        case 'maxiter',
            options.MaxIter = varargin{i+1};
            i = i+1;
        case 'outfile',
            options.outfile = varargin{i+1};
            i = i+1;
        case 'append',
            options.append = true;
        case 'rest_energy',
            options.rest_energy = varargin{i+1};
            i = i+1;
        case 'ebreak',
            options.Ebreak = varargin{i+1};
            i = i+1;
        case 'e0',
            options.E0 = varargin{i+1};
            i = i+1;
        case 'g=gde',
            options.dE_mode = 0;
        case 'trapz',
            options.dE_mode = 1;
        case 'plateau',
            options.dE_mode = 2;
        case 'num_bases',
            options.num_bases = varargin{i+1};
            i = i+1;
    end
    i = i+1;
end

if isempty(options.outfile),
    options.verbose = 0; % silent
else
    options.verbose = 4+options.append; % 4 replace, 5 append
end


function fit = ana_spec_inv(y,dy,Egrid,H0,dt,b,Eout,varargin)

[NT,NY,NE] = spec_inv_check_args(y,dy,Egrid,H0,dt,b);
options = parse_options(varargin{:});
NEout = length(Eout);


if ~options.fxn_bitmap,
    error('No analytical spectrum function(s) selected');
end

% create function table
fxns.pl = struct('bit',options.ASI_FXN_PL,'Nq',2);
fxns.pl.flux = @(q,E)exp(q(1)-q(2)*log(E));
fxns.exp = struct('bit',options.ASI_FXN_EXP,'Nq',2);
fxns.exp.flux = @(q,E)exp(q(1)+q(2)*E);
fxns.rm = struct('bit',options.ASI_FXN_RM,'Nq',2);
fxns.rm.flux = @(q,E)E.*(1+E./options.rest_energy/2).*exp(q(1)+q(2)*E);
fxns.rm2 = struct('bit',options.ASI_FXN_RM2,'Nq',4);
fxns.rm2.flux = @(q,E)E.*(1+E./options.rest_energy/2).*(exp(q(1)+q(2)*E)+exp(q(3)+q(4)*E));
fxns.ple = struct('bit',options.ASI_FXN_PLE,'Nq',2);
fxns.ple.flux = @(q,E,Ebreak,E0)flux_ple;

int_params = int32(zeros(10,1));
real_params = nan(10,1);

int_params(1+0) = NY;
int_params(1+1) = NE;
int_params(1+2) = NEout;
int_params(1+3) = options.fxn_bitmap; % analtyical functions bitmap*/
int_params(1+4) = options.minimizer; % minimizer, 0=BFGS, 3=NM */
int_params(1+5) = options.MaxIter; % maximumn # of iterations */
int_params(1+6) = options.verbose; % 0 = verbose off, no text output, 4 = write messages to outfile */
int_params(1+7) = options.dE_mode;
real_params(1+0) = options.rest_energy;
real_params(1+1) = options.Ebreak;
real_params(1+2) = options.E0;

if options.verbose,
    outFilePtr = libpointer('cstring',options.outfile);
else
    outFilePtr = libpointer('cstring'); % NULL
end
fluxPtr = libpointer('doublePtr',nan(NEout,NT));
lambdaPtr = libpointer('doublePtr',nan(NY,NT));
dlogfluxPtr = libpointer('doublePtr',nan(NEout,NT));
supportPtr = libpointer('doublePtr',nan((options.ASI_MAX_POW2+1)*(2+options.ASI_MAX_NQ+NY),NT));
resultsPtr = libpointer('longPtr',zeros(NT,1));

fit.result = calllib('invlib','ana_spec_inv_multi',NT,y',dy,Egrid,H0',dt,b',int_params,real_params,outFilePtr,Eout,fluxPtr,dlogfluxPtr,lambdaPtr,supportPtr,resultsPtr);

fit.flux = fluxPtr.value';
fit.dlogflux = dlogfluxPtr.value';
fit.lambda = lambdaPtr.value';
fit.result_codes = resultsPtr.value;
support_data = reshape(supportPtr.value,2+options.ASI_MAX_NQ+NY,options.ASI_MAX_POW2+1,NT);
support_data = permute(support_data,[3,2,1]);

functions = fieldnames(fxns);

for k = 1:length(functions),
    f = functions{k};
    if bitand(options.fxn_bitmap,fxns.(f).bit),
        fit.(f) = fxns.(f);
        fit.(f).ell = reshape(support_data(:,k,1),NT,1);
        fit.(f).weight = reshape(support_data(:,k,2),NT,1);
        fit.(f).q = reshape(support_data(:,k,1:fxns.(f).Nq),NT,fxns.(f).Nq);
        fit.(f).lambda = reshape(support_data(:,k,2+options.ASI_MAX_NQ+(1:NY)),NT,NY);
    end
end


function fit = pc_spec_inv(y,dy,Egrid,H0,dt,b,mean_log_flux,basis_vectors,basis_variance,varargin)
% pc_spec_inv: perform principal component spectral inverse
% fit = invlib('pc_spec_inv',y,dy,Egrid,G,dt,b,mean_log_flux,basis_vectors,basis_variance...)
% common arguments have same meaning as for ana_spec_inv
% mean_log_flux - NE x 1 vector of mean log flux
% basis_vectors - NE x Nq matrix of basis vectors (columns)
% basis_variance - Nq x 1 vector of variance in log flux explained by each
%    basis vector
% options: supports appropriate options for ana_spec_inv, plus
% 'num_bases',num_bases - select a reduced number of bases to use

[NT,NY,NE] = spec_inv_check_args(y,dy,Egrid,H0,dt,b);
options = parse_options(varargin{:});

Nq = length(basis_variance);
options.num_bases = min(options.num_bases,Nq);

if numel(mean_log_flux) ~= NE,
    error('mean_log_flux must be NE x 1');
end

if numel(basis_vectors) ~= NE*Nq,
    error('basis_vectors must be NE x Nq');
end


int_params = int32(zeros(10,1));
real_params = nan(10,1);

int_params(1+0) = NY;
int_params(1+1) = NE;
int_params(1+2) = Nq; % size of basis_vectors
int_params(1+3) = options.num_bases; % num bases to actually use
int_params(1+4) = options.minimizer;
int_params(1+5) = options.MaxIter;
int_params(1+6) = options.verbose;
int_params(1+7) = options.dE_mode;

if options.verbose,
    outFilePtr = libpointer('cstring',options.outfile);
else
    outFilePtr = libpointer('cstring'); % NULL
end
fluxPtr = libpointer('doublePtr',nan(NE,NT));
lambdaPtr = libpointer('doublePtr',nan(NY,NT));
dlogfluxPtr = libpointer('doublePtr',nan(NE,NT));
supportPtr = libpointer('doublePtr',nan(1+options.num_bases,NT));
resultsPtr = libpointer('longPtr',zeros(NT,1));

assert(NT==1); % for now, sin
fit.result = calllib('invlib','pc_spec_inv_multi',NT,y',dy,Egrid,H0',dt,b',mean_log_flux,basis_vectors',basis_variance,...
    int_params,real_params,outFilePtr,fluxPtr,dlogfluxPtr,lambdaPtr,supportPtr,resultsPtr);

fit.flux = fluxPtr.value';
fit.dlogflux = dlogfluxPtr.value';
fit.lambda = lambdaPtr.value';
%fit.result_codes = resultsPtr.value;
support_data = supportPtr.value';
fit.ell = support_data(:,1);
fit.q = support_data(:,2:end);


function flux = flux_ple(q,E,Ebreak,E0)
flux = nan(size(E));
f = (E<=Ebreak);
if any(f),
    flux(f) = exp(q(1)-q(2)*E(f));
end
if any(~f),
    flux(~f) = exp(q(1)-q(2)*log(Ebreak)-(E(~f)-E0));
end

function [uniflux,dloguniflux,result] = omni2uni(omniflux,dlogomniflux,method,varargin)

switch(lower(method)),
    case {'tem1'}, method = -1;
    case {'vampola'}, method = -2;
    otherwise
        error('Unknown method "%s"',method);
end

% translate varargin

options.keV = nan;
options.Lm = nan;
options.BB0 = nan;
options.NA = 50;
options.MaxIter = 1e3;
options.root_finder = 0;

i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case {'kev'}, 
            i = i+1;
            options.keV = varargin{i};
        case {'mev'}, 
            i = i+1;
            options.keV = varargin{i}*1e3;
        case {'lm'}, 
            i = i+1;
            options.Lm = varargin{i};
        case {'b/b0'}, 
            i = i+1;
            options.BB0 = varargin{i};
        case {'na'}, 
            i = i+1;
            options.NA = varargin{i};
        case {'maxiter'}, 
            i = i+1;
            options.MaxIter = varargin{i};
        case {'phs'},
            options.root_finder = 0;
        case {'ph'},
            options.root_finder = 1;
        case {'newton'},
            options.root_finder = 2;
        case {'gn'},
            options.root_finder = 3;
        case {'phnhs'},
            options.root_finder = 4;
        case {'phnh'},
            options.root_finder = 5;
        case {'dn'},
            options.root_finder = 6;
        case {'broyden'},
            options.root_finder = 6;
    end
    i = i+1;
end

all_int_params = cell(5,1);
all_real_params = cell(3,1);

all_int_params{1+0} = options.NA;
all_int_params{1+1} = method;
all_int_params{1+2} = 0; % verbose setting = no text output
all_int_params{1+3} = options.root_finder;
all_int_params{1+4} = 1000; % maximumn # of iterations */

all_real_params{1+0} = options.keV;
all_real_params{1+1} = options.BB0;
all_real_params{1+2} = options.Lm;

load_invlib;

int_params = nan(size(all_int_params));
real_params = nan(size(all_real_params));
nullPtr = libpointer('cstring'); % empty pointer to char *
unifluxPtr = libpointer('doublePtr',NaN); % pointer to one double
dlogunifluxPtr = libpointer('doublePtr',NaN); % pointer to one double

N = numel(omniflux);
if numel(dlogomniflux)==1,
    dlogomniflux = repmat(dlogomniflux,size(omniflux));
end
uniflux = nan(size(omniflux));
dloguniflux = nan(size(omniflux));
result = nan(size(omniflux));

for i = 1:N,
    for j = 1:length(int_params),
        if numel(all_int_params{j})==1,
            int_params(j) = all_int_params{j}(1);
        else
            int_params(j) = all_int_params{j}(i);
        end
    end
    for j = 1:length(real_params),
        if numel(all_real_params{j})==1,
            real_params(j) = all_real_params{j}(1);
        else
            real_params(j) = all_real_params{j}(i);
        end
    end
    unifluxPtr.value=NaN;
    dlogunifluxPtr.value=NaN;
    result(i) = calllib('invlib','omni2uni',omniflux(i),dlogomniflux,int_params,real_params,nullPtr,unifluxPtr,dlogunifluxPtr);
    uniflux(i) = unifluxPtr.value;
    dloguniflux(i) = dlogunifluxPtr.value;
end

function load_invlib
% smartly load invlib DLL (.dll or .so)
if ~libisloaded('invlib'),
    if ispc,
        libfile = which('invlib.dll');
    else
        libfile = which('invlib.so');
    end
    hfile = which('invlib.h');
    loadlibrary(libfile,hfile,'alias','invlib');
end

