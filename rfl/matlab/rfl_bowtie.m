function results = rfl_bowtie(inst_info,type,species,exponents,channels,varargin)
% results = rfl_bowtie(inst_info,type,species,exponents,channels)
% results = rfl_bowtie(inst_info,type,species,exponents,channels)
% perform bowtie analysis
% inst_info  - string giving file name or structure containing inst_info
% type - what type of channel to approximate with the bowtie
%   'int' for integral channel
%   'diff' for differential channel
%   'wide' for wide differential channel
%       a wide channel requires an E0 (lower threshold) provided by the
%       'E0' option
% species - string giving species to study (e.g., 'PROT')
% exponents - power law exponents for differential energy spectrum
%   e.g., [3,4,5] for E^-3, E^-4, E^-5
%   exp(-E/T) can be specified with 'T' option (see below).
%   Both power law and exponential types are allowed
%   To do just exponentials, specify 'T' option and leave exponents empty
% channels - cell array of strings, channel names
%  (default or [] does all in inst_info.CHANNEL_NAMES)
% results - struct with members given by channels or inst_info.CHANNEL_NAMES
%   each channel has E0, E1 and G0
%   for diff type, E0 is the effective energy for the channel and
%     E1 is E0
%     G0 is the geometric-efficiency factor times the energy bandwidth (e.g., cm^2 sr MeV)
%   for int type, E0 is the effective energy threshold and
%     E1 is infinity
%     G0 the geometric-efficiency factor (e.g., cm^2 sr)
%   for wide type, E0 is the lower energy threshold
%     E1 is the effective upper energy limit for the channel and
%     G0 is the geometric-efficiency factor (e.g., cm^2 sr)
%   other fields:
%   type - copy of input type
%   fig - figure handle (only present if plotting requested)
% errors on E0 and G0 are given as E0err and G0err
%   these errors are the standard deviation of the log of the E0, G0
%   intersections determined in the bow tie analysis (no sqrt(N) factor)
%
% options:
% 'method','mindG' - Choose E0 or E1 that minimizes spread of G0 
%     Recommended method for determining E0 or E1
%    (default is 'intersection', which is probably not as good)
%     E0err or E1err will be 0
% 'E0',[...] array of same length as CHANNEL_NAMES giving E0 in E_UNIT used
% 'E0','50%' - determine E0 as the first time the response exceeds 50% max
% 'E0','10%' - determine E0 as the first time the response exceeds 10% max
% 'E0','struct' - from the inst_info structure (e.g., inst_info.Elec1.ELE.E0)
%   (if E0 is omitted, 'struct' is assumed, and then '50%' is used if E0 is
%   not in the inst_info structure)
% 'T',T - include exp(-E/T) type exponential spectra for one or more T's.
% 'tail' - struct('E',E_MeV,'fun',@tail)
%    splice high energy tail: j(E) ~tail(E) for E>=E_MeV
%    NOTE: this only works for diff or int bowties
%    and the @tail spectrum is differential
% 'plot' - make diagnostic plot

METHOD_INTERSECT = 1;
METHOD_MINDG = 2;
methods = {'intersect','mindg'};
method = METHOD_INTERSECT;

tail = false;

TYPE_DIFF = 1;
TYPE_WIDE = 2;
TYPE_INT = 3;
types = {'diff','wide','int'};
if (nargin < 2) || ~ismember(lower(type),types),
    error('Must specify type as "int", "diff", or "wide"');
end

if nargin < 3,
    error('Must specify species');
end

if nargin < 4,
    error('Must specify exponents or T''s');
end

if nargin < 5,
    channels = {};
end

type = strmatch(lower(type),types);

E0mode = 'struct';
i = 1;
do_plot = false;
Ts = [];
while i <= length(varargin),
    switch(lower(varargin{i})),
        case {'plot'},
            do_plot=true;
        case {'tail'},
            i = i+1;
            tail = varargin{i};
        case {'t'},
            i = i+1;
            Ts = varargin{i};
        case {'e0'},
            i = i+1;
            E0mode = varargin{i};
        case {'method'},
            i = i+1;
            method = varargin{i};
        otherwise
            error('Unknown option "%s"',varargin{i});
    end
    i = i+1;
end

if ischar(method),
    tmp = find(strcmpi(method,methods));
    if isempty(tmp),
        error('Unknown Method: %d',method);
    end
    method = tmp;
elseif ~ismember(method,1:length(methods)),
    error('Unknown Method: %d',method);
end

if ~isequal(tail,false) && (method ~= METHOD_MINDG),
    error('tail option only supported for method mindG');
end

if ~isequal(tail,false) && ~ismember(type,[TYPE_DIFF,TYPE_INT]),
    error('tail option only supported for type int or diff');
end

clear sinfo

sinfo.PL = struct('type','PL','index',1,'linestyle','k-','special',false);
sinfo.EXP = struct('type','EXP','index',2,'linestyle','k--','special',false);

Nn = length(exponents);
NT = length(Ts);
Ns = NT+Nn;
spectra = cell(Ns,1);
for i = 1:Nn,
    spectra{i} = setfield(sinfo.PL,'param',exponents(i));
    n = exponents(i);
    spectra{i}.j = @(E)flux_PL(E,n,tail);
    spectra{i}.special = (n<=1) || (type==TYPE_WIDE);
    switch(type),
        case TYPE_INT,
            spectra{i}.G = @(c,E,E0,Emax)G_PL_int_tail(c,E,E0,Emax,n,tail);
        case TYPE_DIFF,
            %spectra{i}.G = @(c,E,E0,Emax)c*E.^n;
            spectra{i}.G = @(c,E,E0,Emax)G_PL_diff_tail(c,E,E0,Emax,n,tail);
        case TYPE_WIDE,            
            if n==1, % need to include Emax and special case n
                spectra{i}.G = @(c,E,E0,Emax)c./(log(E)-log(E0));
            else % treate Emax as infinity
                spectra{i}.G = @(c,E,E0,Emax)c*(n-1)./(E0^(1-n)-E.^(1-n));
            end
    end
    
end
for i = 1:NT,
    T = Ts(i);
    spectra{Nn+i} = setfield(sinfo.EXP,'param',T);
    spectra{Nn+i}.j = @(E)flux_EXP(E,T,tail);
    switch(type),
        case TYPE_INT,
            %spectra{Nn+i}.G = @(c,E,E0,Emax)c.*exp(E./T)./T;
            spectra{Nn+i}.G = @(c,E,E0,Emax)G_EXP_int_tail(c,E,E0,Emax,T,tail);
        case TYPE_DIFF,
            %spectra{Nn+i}.G = @(c,E,E0,Emax)c.*exp(E./T);
            spectra{Nn+i}.G = @(c,E,E0,Emax)G_EXP_diff_tail(c,E,E0,Emax,T,tail);
        case TYPE_WIDE,
            spectra{Nn+i}.G = @(c,E,E0,Emax)c./T./(exp(-E0./T)-exp(-E./T));
    end
end

% build intersect function table
intersect_func = cell(2,2);
switch(type),
    case TYPE_INT,
        intersect_func{sinfo.PL.index,sinfo.PL.index} = @intersect_INT_PL_PL;
        intersect_func{sinfo.PL.index,sinfo.EXP.index} = @intersect_INT_PL_EXP;
        intersect_func{sinfo.EXP.index,sinfo.PL.index} = @intersect_INT_EXP_PL;
        intersect_func{sinfo.EXP.index,sinfo.EXP.index} = @intersect_INT_EXP_EXP;
    case TYPE_DIFF,
        intersect_func{sinfo.PL.index,sinfo.PL.index} = @intersect_DIFF_PL_PL;
        intersect_func{sinfo.PL.index,sinfo.EXP.index} = @intersect_DIFF_PL_EXP;
        intersect_func{sinfo.EXP.index,sinfo.PL.index} = @intersect_DIFF_EXP_PL;
        intersect_func{sinfo.EXP.index,sinfo.EXP.index} = @intersect_DIFF_EXP_EXP;
    case TYPE_WIDE,
        intersect_func{sinfo.PL.index,sinfo.PL.index} = @intersect_WIDE_PL_PL;
        intersect_func{sinfo.PL.index,sinfo.EXP.index} = @intersect_WIDE_PL_EXP;
        intersect_func{sinfo.EXP.index,sinfo.PL.index} = @intersect_WIDE_EXP_PL;
        intersect_func{sinfo.EXP.index,sinfo.EXP.index} = @intersect_WIDE_EXP_EXP;
end

if ischar(inst_info),
    inst_info = rfl_load_inst_info(inst_info);
end

if isempty(channels),
    channels = inst_info.CHANNEL_NAMES;
end

if isnumeric(E0mode) && (numel(E0mode) ~= numel(channels)),
    error('Supplied E0 array must have one entry per channel (%d)',numel(channels));
end

clear results
for ichan = 1:length(channels),
    chan = channels{ichan};
    if ~isfield(inst_info,chan),
        error('Channel %s not found in inst_info',chan);
    end
    results.(chan) = struct('G0',NaN,'E0',NaN,'E1',NaN,'type',types{type});
    ispec = find(strcmpi(species,inst_info.(chan).SPECIES));
    if isempty(ispec),
        warning('Species %s not found in %s',species,chan);
        continue
    end
    sp = inst_info.(chan).SPECIES{ispec};
    resp = inst_info.(chan).(sp);
    [hE,result_code] = resp.make_hE(resp,resp.E_GRID,[]);
    R = hE./rfl_make_deltas(resp.E_GRID);
    if type == TYPE_WIDE,
        if isnumeric(E0mode),
            E0 = E0mode(ichan);
        else
            if strcmpi(E0mode,'struct'),
                if isfield(resp,'E0'),
                    E0 = resp.E0;
                else
                    fprintf('E0 not supplied for %s, defaulting to 50%% method\n',chan);
                    E0mode = '50%';
                end
            end
            if strcmp(E0mode,'50%'),
                E0 = resp.E_GRID(min(find(R>=max(R)/2)));
            end
            if strcmp(E0mode,'10%'),
                E0 = resp.E_GRID(min(find(R>=max(R)/10)));
            end
        end
        E_GRID = resp.E_GRID(resp.E_GRID>=E0);
    else
        E0 = nan;
        E_GRID = resp.E_GRID;
    end
    Emax = E_GRID(end);
    
    if do_plot,
        fig = figure;
        results.(chan).fig = fig;
    end
    
    C = nan(Ns,1);
    % find reference counts
    for is = 1:Ns,
        C(is) = hE(:)'*spectra{is}.j(resp.E_GRID(:));
    end
    % find intersections
    intersect = [];
    if method == METHOD_INTERSECT,
        for ipass = 1:2,
            % pass 1: do same type
            % pass 2: do hybrid type w/ guess since these can have multiple crossings
            if ipass == 1,
                Eguess = nan;
            else
                Eguess = median(intersect(:,1)); % initial guess from same-type intersections
            end
            for is1 = 1:Ns,
                s1 = spectra{is1};
                c1 = C(is1);
                if do_plot && (ipass==1),
                    figure(fig);
                    loglog(E_GRID,s1.G(c1,E_GRID,E0,Emax),s1.linestyle,'tag',sprintf('%s: %g',spectra{is}.type,spectra{is}.param));
                    hold on;
                end
                for is2 = (is1+1):Ns,
                    s2 = spectra{is2};
                    ishybrid = isequal(s1.type,s2.type);
                    if (ipass==1) == (ishybrid || s1.special || s2.special), % do this on pass 1 for same types, on pass 2 for hybrid types
                        c2 = C(is2);
                        E = intersect_func{s1.index,s2.index}(s1,s2,c1,c2,E0,Emax,Eguess);
                        G0 = s1.G(c1,E,E0,Emax);
                        if (E<=0) || ~isreal(E) || ~isfinite(E) || (G0<=0) || ~isreal(G0) || ~isfinite(G0),
                            warning('No real, positive, finite fit found for %s/%g - %s/%g',s1.type,s1.param,s2.type,s2.param);
                            continue;
                        end
                        intersect = [intersect;E G0 is1 is2];
                    end
                end
            end
        end
    else % method == METHOD_MINDG
        G0 = nan(length(E_GRID),Ns);
        for is = 1:Ns,
            G0(:,is) = spectra{is}.G(C(is),E_GRID,E0,Emax);
            if do_plot,
                figure(fig);
                loglog(E_GRID,G0(:,is),spectra{is}.linestyle,'tag',sprintf('%s: %g',spectra{is}.type,spectra{is}.param));
                hold on;
            end
        end
        [stdG,imin] = min(std(log(G0),[],2)); % find best guess
        imin = imin(1);
        intersect = nan(Ns,4);
        intersect(:,1) = E_GRID(imin);
        intersect(:,2) = G0(imin,:)';
        intersect(:,3) = 1:Ns;
        intersect(:,4) = 1:Ns;
    end
    if do_plot,
        figure(fig);
        loglog(intersect(:,1),intersect(:,2),'bo');
    end
    solution = exp(median(log(intersect(:,1:2))));
    sol_error = std(log(intersect(:,1:2)));
    results.(chan).G0 = solution(2);
    results.(chan).G0err = sol_error(2);
    switch(type),
        case TYPE_INT,
            Evar = 'E0';
            results.(chan).E0 = solution(1);
            results.(chan).E0err = sol_error(1);
            results.(chan).E1 = inf;
            results.(chan).E1err = 0;
        case TYPE_DIFF,
            Evar = 'E0';
            results.(chan).E0 = solution(1);
            results.(chan).E0err = sol_error(1);
            results.(chan).E1 = results.(chan).E0;
            results.(chan).E1err = results.(chan).E0err;
        case TYPE_WIDE,
            Evar = 'E1';
            results.(chan).E1 = solution(1);
            results.(chan).E1err = sol_error(1);
            results.(chan).E0 = E0;
            results.(chan).E0err = 0;
    end
    
    if do_plot,
        figure(fig);
        if any(R(:)>0),
            loglog(solution(1),solution(2),'rx','linew',3);
            axmin = min(intersect(:,1:2));
            axmax = max(intersect(:,1:2));
            axmin(1) = min(axmin(1),E_GRID(1));
            axmax(1) = max(axmax(1),E_GRID(end));
            axmin(2) = min(axmin(2),min(R(R>0)));
            axmax(2) = max(axmax(2),max(R));
            axmin = 10.^(floor(-0.5+log10(axmin)));
            axmax = 10.^(ceil(0.5+log10(axmax)));
            loglog(resp.E_GRID,max(R,axmin(2)),'g-','linew',2); % epsdEG vs E
            axis([axmin(1) axmax(1) axmin(2) axmax(2)]);
        else
            axis([0,1,0,1]);
            text(0.5,0.5,'No positive entries in response','horiz','center','vert','mid');
        end
        switch(type),
            case {TYPE_INT,TYPE_WIDE},
                Gsym = 'G0';
                G_UNIT = sprintf('%s^2sr',resp.L_UNIT);
            case TYPE_DIFF,
                Gsym = 'G0dE';
                G_UNIT = sprintf('%s%s^2sr',resp.E_UNIT,resp.L_UNIT);
        end
        title(sprintf('%s,%s : %s=%g (%.1f%%), %s=%g (%.1f%%)',...
            chan,types{type},Evar,solution(1),sol_error(1)*100,...
            Gsym,results.(chan).G0,results.(chan).G0err*100),'interpreter','none');
        xlabel(sprintf('%s, %s',Evar,resp.E_UNIT));
        ylabel(sprintf('%s, %s',Gsym,G_UNIT));
        grid on;
    end
end

% intersect functions all have same syntax:
% input two spectra and two counts
% and E0 which is often ignored
% Eguess is also often ignored
function E = intersect_INT_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)
n1 = s1.param;
n2 = s2.param;
if n1 > n2, % try again with n1 < n2
    E = intersect_INT_PL_PL(s2,s1,c2,c1,E0,Emax,Eguess);
    return
end
if n1<1, % have to deal with Emax
    % non-analytic cases
    % n2 < 1 cases
    % G0 = c1*(n1-1)/(E^(n1-1)-Emax^(n1-1))
    % G0 = c2*(n2-1)/(E^(n2-1)-Emax^(n2-1))
    % c1*(n1-1)/(E^(n1-1)-Emax(n2-1)) = c2*(n2-1)/(E^(n2-1)-Emax^(n1-1))
    % n2 = 1, n2>1 cases also not analytic
    if ~isfinite(Eguess),
        Eguess = Emax*0.95;
    end
    E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess); % non-analytic
elseif n1 == 1, % n1 = 1, n2 > 1
    
    % G0 = c1/(log(Emax)-log(E))
    % G0 = c2*E0^(n2-1)*(n2-1)
    % c1*E0^(n1-1)*(n1-1) = c1/(log(Emax)-log(E))
    % no analytical solution
    if ~isfinite(Eguess), % guess at n=1.01 or half-way to n2
        n1 = min(1.01,(n1+n2)/2);
        Eguess = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1)); % E0
    end
    E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess); % non-analytic
else
    % G0 = c1*E0^(n1-1)*(n1-1)
    % G0 = c2*E0^(n2-1)*(n2-1)
    % c1*E0^(n1-1)*(n1-1) = c2*E0^(n2-1)*(n2-1)
    % c1/c2*(n1-1)/(n2-1) = E0^(n2-n1)
    % E0 = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1))
    E = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1)); % E0
end

function E = intersect_INT_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
% c2/T2*exp(E/T2) = (n1-1)*c1*E^(n1-1)
E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess);

function E = intersect_INT_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)
E = intersect_INT_PL_EXP(s2,s1,c2,c1,E0,Emax,Eguess);

function E = intersect_INT_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
T1 = s1.param;
T2 = s2.param;
E = (log(c2)-log(T2)-log(c1)+log(T1))/(1/T1-1/T2);

function E = intersect_DIFF_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)
n1 = s1.param;
n2 = s2.param;
% G0 = c1*E0^n1
% G0 = c2*E0^n2
% c1*E0^n1 = c2*E0^n2
% c1/c2 = E0^(n2-n1)
% E0 = (c1/c2)^(1/(n2-n1))
E = (c1/c2)^(1/(n2-n1)); % E0

function E = intersect_DIFF_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
% log(c1)-log(c2) = E/T2-n1*log(E)
E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess);

function E = intersect_ID_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
[logE,fval,exitflag] = fzero(@(logE)log(s1.G(c1,exp(logE),E0,Emax))-log(s2.G(c2,exp(logE),E0,Emax)),log(Eguess));
if exitflag==1,
    E = exp(logE);
else
    E = nan;
end

function E = intersect_DIFF_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)
E = intersect_DIFF_PL_EXP(s2,s1,c2,c1);

function E = intersect_DIFF_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
T1 = s1.param;
T2 = s2.param;
E = (log(c1)-log(c2))/(1/T2-1/T1);

function E = intersect_WIDE_PL_PL(s1,s2,c1,c2,E0,Emax,Eguess)
% G0 = c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)]
% G0 = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
% c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)] = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
% no analytical solution, solve numerically
E = intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess);

function E = intersect_WIDE_PL_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
% no analytical solution
E = intersect_WIDE(s2,s1,c2,c1,E0,Emax,Eguess);

function E = intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess) % generic
if isfinite(Eguess),
    logdE = log(Eguess-E0);
else
    logdE = log(E0)-2;
end
[logdE,fval,exitflag] = fzero(@(logdE)log(s1.G(c1,E0+exp(logdE),E0,Emax))-log(s2.G(c2,E0+exp(logdE),E0,Emax)),logdE);
if exitflag==1,
    E = E0+exp(logdE);
else
    E = nan;
end

function E = intersect_WIDE_EXP_PL(s1,s2,c1,c2,E0,Emax,Eguess)
E = intersect_WIDE_PL_EXP(s2,s1,c2,c1,E0,Emax,Eguess);

function E = intersect_WIDE_EXP_EXP(s1,s2,c1,c2,E0,Emax,Eguess)
% c1*(E-E0)/T1/(exp(-E0/T1)-exp(-E/T1)) = c2*(E-E0)/T2/(exp(-E0/T2)-exp(-E/T2));
E = intersect_WIDE(s1,s2,c1,c2,E0,Emax,Eguess);


function G = G_PL_diff_tail(c,E,E0,Emax,n,tail)
G = c*E.^n;
if ~isequal(tail,false),
    iE = E>tail.E;
    if any(iE),
        G(iE) = c./tail.fun(E(iE))*tail.E^n*tail.fun(tail.E);
    end
end

function G = G_EXP_diff_tail(c,E,E0,Emax,T,tail)
G = c.*exp(E./T);
if ~isequal(tail,false),
    iE = E>tail.E;
    if any(iE),
        G(iE) = c./tail.fun(E(iE))*exp(tail.E/T)*tail.fun(tail.E);
    end
end

function G = G_PL_int_tail(c,E,E0,Emax,n,tail)
% jdiff = E^-n
% G = c/Jint
Elim = Emax;
if ~isequal(tail,false),
    if tail.E>Emax,
        tail = false;
    else
        Elim = tail.E; % Elim serves as Emax
    end
end
% integrals up to Elim~Emax or tail.E
if n<1, 
    %G = c*(n-1)./(Emax.^(n-1)-E.^(n-1));
    Jint = (Elim.^(n-1)-E.^(n-1))/(n-1);
elseif n==1, 
    %G = c./(log(Emax)-log(E));
    Jint = log(Elim)-log(E);
else 
    % G = c*E.^(n-1)*(n-1);
    Jint = (E.^(1-n)-Elim.^(1-n))/(n-1);
end
% now do tail integral if needed
if ~isequal(tail,false),
    J0int = Jint;
    jtail = tail.fun(tail.E);
    j0Etail = tail.E^-n;
    itail = E>tail.E;
    if any(~itail),
        JprimeEtail = j0Etail/jtail*quad(tail.fun,tail.E,Emax);
        Jint(~itail) = J0int(~itail)+JprimeEtail;
    end
    if any(itail),
        for i = find(itail(:))'
            Jint(i) = j0Etail/jtail*quad(tail.fun,E(i),Emax);
        end
    end
end

G = c./Jint;

function G = G_EXP_int_tail(c,E,E0,Emax,T,tail)
% jdiff = exp(-E/T)
% G = c/Jint
% G = c.*exp(E./T)./T;
Jint = T.*exp(-E./T);
% now do tail integral if needed
if ~isequal(tail,false),
    J0int = Jint;
    jtail = tail.fun(tail.E);
    j0Etail = exp(-tail.E./T);
    itail = E>tail.E;
    if any(~itail),
        JprimeEtail = j0Etail/jtail*quad(tail.fun,tail.E,Emax);
        Jint(~itail) = J0int(~itail)+JprimeEtail;
    end
    if any(itail),
        for i = find(itail(:))'
            Jint(i) = j0Etail/jtail*quad(tail.fun,E(i),Emax);
        end
    end
end

G = c./Jint;

function flux = flux_PL(E,n,tail)
flux = E.^(-n);
if ~isequal(tail,false),
    itail = E>tail.E;
    if any(itail),
        flux(itail) = tail.fun(E(itail))/tail.fun(tail.E)*tail.E.^(-n);
    end
end

function flux = flux_EXP(E,T,tail)
flux = exp(-E./T);
if ~isequal(tail,false),
    itail = E>tail.E;
    if any(itail),
        flux(itail) = tail.fun(E(itail))/tail.fun(tail.E)*exp(-tail.E./T);
    end
end
