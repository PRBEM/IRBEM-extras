function varargout = invlib(what,varargin)
% varargout = invlib(what,varargin)
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
% set rest_energy for RM and RM2: ...,'rest_energy',rest_energy,...
% set Ebreak for PLE: ...,'Ebreak',Ebreak,...
% set E0 for PLE: ...,'E0',E0,...
% set maximum fit iterations: ...,'MaxIter',MaxIter,...
% set verbose output text file: ...,'outfile',filename,...
% set append to output file: ...,'append',...
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

if ~libisloaded('invlib'),
    if ispc,
        libfile = which('invlib.dll');
    else
        libfile = which('invlib.so');
    end
    hfile = which('invlib.h');
    loadlibrary(libfile,hfile,'alias','invlib');
end

switch(lower(what)),
    case {'ana_spec_inv'},
        [varargout{:}] = ana_spec_inv(varargin{:});
    otherwise
        error('%s: Unknown "what": "%s"',mfilename,what);
end

function fit = ana_spec_inv(y,dy,Egrid,H0,dt,b,Eout,varargin)

[NT,NY] = size(y);
NE = length(Egrid);
NEout = length(Eout);
if numel(dy) ~= dy,
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

% constants from specinv.h
ASI_FXN_PL   = (1);
ASI_FXN_EXP  = (2);
ASI_FXN_RM   = (4);
ASI_FXN_PLE  = (8);
ASI_FXN_RM2  = (16);
ASI_MAX_POW2 = 4;
ASI_MAX_NQ = 10;

% constants from optim.h
OPTIM_MIN_BFGS  = (0);
OPTIM_MIN_FR    = (1);
OPTIM_MIN_PR    = (2);
OPTIM_MIN_NM    = (3);

minimizer = 0; % BFGS
MaxIter = 1000;
fxn_bitmap = 0;
outfile = '';
rest_energy = 0.511; % electron rest energy, MeV
Ebreak = 100; % power-law to exponential tail transition energy, MeV
E0 = 345; % power-law to exponential tail transition energy, MeV
append = false; % append to text file?

i = 1;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case 'bfgs',
            minimizer = OPTIM_MIN_BFGS;
        case 'fr',
            minimizer = OPTIM_MIN_FR;
        case 'pr',
            minimizer = OPTIM_MIN_PR;
        case 'nm',
            minimizer = OPTIM_MIN_NM;
        case 'pl',
            fxn_bitmap = bitor(fxn_bitmap,ASI_FXN_PL);
        case 'exp',
            fxn_bitmap = bitor(fxn_bitmap,ASI_FXN_EXP);
        case 'rm',
            fxn_bitmap = bitor(fxn_bitmap,ASI_FXN_RM);
        case 'ple',
            fxn_bitmap = bitor(fxn_bitmap,ASI_FXN_PLE);
        case 'rm2',
            fxn_bitmap = bitor(fxn_bitmap,ASI_FXN_RM2);
        case 'maxiter',
            MaxIter = varargin{i+1};
            i = i+1;
        case 'outfile',
            outfile = varargin{i+1};
            i = i+1;
        case 'append',
            append = true;
        case 'rest_energy',
            rest_energy = varargin{i+1};
            i = i+1;
        case 'ebreak',
            Ebreak = varargin{i+1};
            i = i+1;
        case 'e0',
            E0 = varargin{i+1};
            i = i+1;
    end
    i = i+1;
end

if ~fxn_bitmap,
    error('No analytical spectrum function(s) selected');
end

if isempty(outfile),
    verbose = 0; % silent
else
    verbose = 4+append; % 4 replace, 5 append
end

% create function table
fxns.pl = struct('bit',ASI_FXN_PL,'Nq',2);
fxns.pl.flux = @(q,E)exp(q(1)-q(2)*log(E));
fxns.exp = struct('bit',ASI_FXN_EXP,'Nq',2);
fxns.exp.flux = @(q,E)exp(q(1)+q(2)*E);
fxns.rm = struct('bit',ASI_FXN_RM,'Nq',2);
fxns.rm.flux = @(q,E)E.*(1+E./rest_energy/2).*exp(q(1)+q(2)*E);
fxns.rm2 = struct('bit',ASI_FXN_RM2,'Nq',4);
fxns.rm2.flux = @(q,E)E.*(1+E./rest_energy/2).*(exp(q(1)+q(2)*E)+exp(q(3)+q(4)*E));
fxns.ple = struct('bit',ASI_FXN_PLE,'Nq',2);
fxns.ple.flux = @(q,E,Ebreak,E0)flux_ple;

int_params = int32(zeros(7,1));
real_params = nan(3,1);

int_params(1+0) = NC;
int_params(1+1) = NE;
int_params(1+2) = NEout;
int_params(1+3) = fxn_bitmap; % analtyical functions bitmap*/
int_params(1+4) = minimizer; % minimizer, 0=BFGS, 3=NM */
int_params(1+5) = MaxIter; % maximumn # of iterations */
int_params(1+6) = verbose; % 0 = verbose off, no text output, 4 = write messages to outfile */
real_params(1+0) = rest_energy;
real_params(1+1) = Ebreak;
real_params(1+2) = E0;

outFilePtr = libpointer('cstring',outfile);
fluxPtr = libpointer('doublePtr',nan(NEout,NT));
lambdaPtr = libpointer('doublePtr',nan(NY,NT));
dlogfluxPtr = libpointer('doublePtr',nan(NEout,NT));
supportPtr = libpointer('doublePtr',nan((ASI_MAX_POW2+1)*(2+ASI_MAX_NQ+NY)),NT);
resultsPtr = libpointer('int32Ptr',zeros(NT,1));

% set up weights for trapezoidal energy integral
dE = [(Egrid(2)-Egrid(1)) ; Egrid(3:end)-Egrid(1:(end-2)) ; Egrid(end)-Egrid(end-1)]/2;
dE(1) = [dE(1)/2;dE(2:end);dE(end)/2];

H0 = G*spdiag(dE,0,NE,NE);

fit.result = calllib('invlib','ana_spec_inv_multi',y',dy,Egrid,H0',dt,b',int_params,real_params,outFilePtr,Eout,fluxPtr,dlogfluxPtr,lambdaPtr,supportPtr,resultsPtr);

fit.flux = fluxPtr.value';
fit.dlogflux = dlogfluxPtr.value';
fit.lambda = lambdaPtr.value';
fit.result_codes = resultsPtr.value;
support_data = reshape(supportPtr.value,2+ASI_MAX_NQ+NY,ASI_MAX_POW2+1,NT);
support_data = permute(support_data,[3,2,1]);

functions = fieldnames(fxns);

for k = 1:length(functions),
    f = functions{k};
    if bitand(fxn_bitmap,fxns.(f).bit),
        fit.(f) = fxns.(f);
        fit.(f).ell = support_data(:,k,1);
        fit.(f).weight = support_data(:,k,2);
        fit.(f).q = support_data(:,k,1:fxns.(f).Nq);
        fit.(f).lambda = support_data(:,k,2+ASI_MAX_NQ+(1:NY));
    end
end

function flux = flux_ple(q,E,Ebreak,E0)
flux = nan(size(E));
f = (E<=Ebreak);
if any(f),
    flux(f) = exp(q(1)-q(2)*E(f));
end
if any(~f),
    flux(~f) = exp(q(1)-q(2)*log(Ebreak)-(E(~f)-E0));
end
