function results = rfl_bowtie(inst_info,type,species,exponents,channels,varargin)
% results = rfl_bowtie(inst_info,type,species,exponents,channels)
% perform bowtie analysis
% inst_info  - string giving file name or structure containing inst_info
% type - what type of channel to approximate with the bowtie
%   'int' for integral channel
%   'diff' for differential channel
%   'wide' for wide differential channel
%       a wide channel requires an E0 (lower threshold) provided a priori
%       in an optional input 'E0' or 
%       in the inst_info structure for each channel
% species - string giving species to study (e.g., 'PROT')
% exponents - power law exponents for differential energy spectrum
%   e.g., [3,4,5] for E^-3, E^-4, E^-5
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
%   for wide type, E0 is copied from the input inst_info
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
% 'E0',E0 - array of same length as CHANNEL_NAMES giving E0 in E_UNIT used
%   by the inst_info structure
% 'plot' - make diagnostic plot

if ischar(inst_info),
    inst_info = rfl_load_inst_info(inst_info);
end

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
    error('Must specify exponents');
end

if nargin < 5,
    channels = {};
end

type = strmatch(lower(type),types);

E0s = [];
i = 1;
do_plot = false;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case {'plot'},
            do_plot=true;
        case {'e0'},
            i = i+1;
            E0s = varargin{i};
        otherwise
            error('Unknown option "%s"',varargin{i});
    end
    i = i+1;
end

if isempty(channels),
    channels = inst_info.CHANNEL_NAMES;
end

if ~isempty(E0s) && (numel(E0s) ~= numel(channels)),
    error('Supplied E0 array must have one entry per channel (%d)',numel(channels));
end

Nn = length(exponents);
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
    if type == TYPE_WIDE,
        if ~isempty(E0s),
            E0 = E0s(ichan);
        else
            if isfield(resp,'E0'),
                E0 = resp.E0;
            else
                error('E0 missing from %s (required for "wide" fit)',chan);
            end
        end
        E_GRID = resp.E_GRID(resp.E_GRID>=E0);
    else
        E_GRID = resp.E_GRID;
    end
    [hE,result_code] = resp.make_hE(resp,resp.E_GRID,[]);
    
    if do_plot,
        fig = figure;
        results.(chan).fig = fig;
    end
    
    C = nan(Nn,1);
    % find reference counts
    for iexp = 1:Nn,
        n = exponents(iexp);
        j = resp.E_GRID(:).^-n;
        c = hE(:)'*j;
        % for diff:
        %   c = E0^-n*G0
        %   G0 = c*E0^n
        % for int:
        %   c = E0^(1-n)/(n-1)*G0
        %   G0 = c*E0^(n-1)*(n-1)
        % for wide:
        %   c = [E0^(1-n)-E1^(1-n)]/(n-1)*G0
        %   G0 = c*(n-1)/[E0^(1-n)-E1^(1-n)]
        C(iexp) = c;
    end
    % find intersections
    intersect = [];
    for iexp1 = 1:Nn,
        n1 = exponents(iexp1);
        c1 = C(iexp1);
        if do_plot,
            figure(fig);
            switch(type),
                case TYPE_INT,
                    loglog(E_GRID,c1*E_GRID.^(n1-1)*(n1-1),'k-'); % G0 vs E0
                case TYPE_DIFF,
                    loglog(E_GRID,c1*E_GRID.^n1,'k-'); % G0 vs E0
                case TYPE_WIDE,
                    loglog(E_GRID,c1*(n1-1)./(E0^(1-n1)-E_GRID.^(1-n1)),'k-'); % G0 vs E0
            end
            hold on;
        end
        for iexp2 = (iexp1+1):Nn,
            n2 = exponents(iexp2);
            c2 = C(iexp2);
            switch(type),
                case TYPE_INT,
                    % G0 = c1*E0^(n1-1)*(n1-1)
                    % G0 = c2*E0^(n2-1)*(n2-1)
                    % c1*E0^(n1-1)*(n1-1) = c2*E0^(n2-1)*(n2-1)
                    % c1/c2*(n1-1)/(n2-1) = E0^(n2-n1)
                    % E0 = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1))
                    E = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1)); % E0
                    G0 = c1*E^(n1-1)*(n1-1);
                case TYPE_DIFF,
                    % G0 = c1*E0^n1
                    % G0 = c2*E0^n2
                    % c1*E0^n1 = c2*E0^n2
                    % c1/c2 = E0^(n2-n1)
                    % E0 = (c1/c2)^(1/(n2-n1))
                    E = (c1/c2)^(1/(n2-n1)); % E0
                    G0 = c1*E^n1;
                case TYPE_WIDE,
                    % G0 = c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)]
                    % G0 = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
                    % c1*(n1-1)/[E0^(1-n1)-E1^(1-n1)] = c2*(n2-1)/[E0^(1-n2)-E1^(1-n2)]
                    % no analytical solution, solve numerically
                    G0fun = @(E1,c,n)c*(n-1)./(E0.^(1-n)-E1.^(1-n));
                    %E =
                    %fzero(@(E1)G0fun(E1,c1,n1)-G0fun(E1,c2,n2),1.5*E0); % poor convergence
                    %E =
                    %fminsearch(@(E1)(G0fun(E1,c1,n1)-G0fun(E1,c2,n2)).^2,1.5*E0); % poor convergence
                    % instead define E1 = E0+exp(dE1) and search in dE1
                    dE1 = fzero(@(dE1)log(G0fun(E0+exp(dE1),c1,n1))-log(G0fun(E0+exp(dE1),c2,n2)),-2);
                    E = E0+exp(dE1);
                    G0 = c1*(n1-1)/(E0^(1-n1)-E^(1-n1));
            end
            if do_plot,
                figure(fig);
                loglog(E,G0,'bo');
            end
            intersect = [intersect;E G0];
        end
    end
    solution = exp(median(log(intersect)));
    sol_error = std(log(intersect));
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
        loglog(solution(1),solution(2),'rx','linew',3);
        yeff = hE./rfl_make_deltas(resp.E_GRID);
        axmin = min(intersect);
        axmax = max(intersect);
        axmin(1) = min(axmin(1),E_GRID(1));
        axmax(1) = max(axmax(1),E_GRID(end));
        axmin(2) = min(axmin(2),min(yeff(yeff>0)));
        axmax(2) = max(axmax(2),max(yeff));
        axmin = 10.^(floor(-0.5+log10(axmin)));
        axmax = 10.^(ceil(0.5+log10(axmax)));
        loglog(resp.E_GRID,max(yeff,axmin(2)),'g-','linew',2); % epsdEG vs E
        axis([axmin(1) axmax(1) axmin(2) axmax(2)]);
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
