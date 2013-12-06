function results = rfl_bowtie(inst_info,type,species,exponents,channels,varargin)
% results = rfl_bowtie(inst_info,type,species,exponents,channels)
% perform bowtie analysis
% inst_info  - string giving file name or structure containing inst_info
% type - 'int' for integral channel or 'diff' for differential channel
% species - string giving species to study (e.g., 'PROT')
% exponents - power law exponents for differential energy spectrum
%   e.g., [3,4,5] for E^-3, E^-4, E^-5
% channels - cell array of strings, channel names
%  (default or [] does all in inst_info.CHANNEL_NAMES)
% results - struct with members given by channels or inst_info.CHANNEL_NAMES
%   each channel has E0 and G0 - the nominal energy E0 and
%   for diff type, G0 includes the energy bandwidth (e.g., cm^2 sr MeV)
%   for int type, G0 includes only the geometric factor (e.g., cm^2 sr)
%   other fields: 
%   type - copy of input type
%   fig - figure handle (only present if plotting requested)
% errors on E0 and G0 are given as E0err and G0err
%   these errors are the standard deviation of the log of the E0, G0
%   intersections determined in the bow tie analysis (no sqrt(N) factor)
%
% options:
% 'plot' - make diagnostic plot

if ischar(inst_info),
    inst_info = rfl_load_inst_info(inst_info);
end

if (nargin < 2) || ~ismember(lower(type),{'int','diff'}),
    error('Must specify type as "int" or "diff"');
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

i = 1;
do_plot = false;
while i <= length(varargin),
    switch(lower(varargin{i})),
        case {'plot'},
            do_plot=true;
        otherwise
            error('Unknown option "%s"',varargin{i});
    end
    i = i+1;
end

if isempty(channels),
    channels = inst_info.CHANNEL_NAMES;
end

intE = strcmpi(type,'int');

Nn = length(exponents);
clear results
for ichan = 1:length(channels),
    chan = channels{ichan};
    if ~isfield(inst_info,chan),
        error('Channel %s not found in inst_info',chan);
    end
    results.(chan) = struct('G0',NaN,'E0',NaN,'type',lower(type));
    ispec = find(strcmpi(species,inst_info.(chan).SPECIES));
    if isempty(ispec),
        warning('Species %s not found in %s',species,chan);
        continue
    end
    sp = inst_info.(chan).SPECIES{ispec};
    resp = inst_info.(chan).(sp);
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
        C(iexp) = c;
    end
    % find intersections
    intersect = [];
    for iexp1 = 1:Nn,
        n1 = exponents(iexp1);
        c1 = C(iexp1);
        if do_plot,
            figure(fig);
            if intE,
                loglog(resp.E_GRID,c1*resp.E_GRID.^(n1-1)*(n1-1),'k-'); % G0 vs E0
            else
                loglog(resp.E_GRID,c1*resp.E_GRID.^n1,'k-'); % G0 vs E0
            end
            hold on;
        end
        for iexp2 = (iexp1+1):Nn,
            n2 = exponents(iexp2);
            c2 = C(iexp2);
            if intE,
                % G0 = c1*E0^(n1-1)*(n1-1)
                % G0 = c2*E0^(n2-1)*(n2-1)
                % c1*E0^(n1-1)*(n1-1) = c2*E0^(n2-1)*(n2-1)
                % c1/c2*(n1-1)/(n2-1) = E0^(n2-n1)
                % E0 = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1))
                E0 = (c1/c2*(n1-1)/(n2-1))^(1/(n2-n1));
                G0 = c1*E0^(n1-1)*(n1-1);
            else
                % G0 = c1*E0^n1
                % G0 = c2*E0^n2
                % c1*E0^n1 = c2*E0^n2
                % c1/c2 = E0^(n2-n1)
                % E0 = (c1/c2)^(1/(n2-n1))
                E0 = (c1/c2)^(1/(n2-n1));
                G0 = c1*E0^n1;
            end
            if do_plot,
                figure(fig);                
                loglog(E0,G0,'bo');
            end
            intersect = [intersect;E0 G0];
        end
    end
    solution = exp(median(log(intersect)));
    sol_error = std(log(intersect));
    results.(chan).E0 = solution(1);
    results.(chan).G0 = solution(2);
    results.(chan).E0err = sol_error(1);
    results.(chan).G0err = sol_error(2);
    if do_plot,
        figure(fig);
        loglog(solution(1),solution(2),'rx','linew',3);
        yeff = hE./rfl_make_deltas(resp.E_GRID);
        axmin = min(intersect);
        axmax = max(intersect);
        axmin(1) = min(axmin(1),resp.E_GRID(1));
        axmax(1) = max(axmax(1),resp.E_GRID(end));
        axmin(2) = min(axmin(2),min(yeff(yeff>0)));
        axmax(2) = max(axmax(2),max(yeff));
        axmin = 10.^(floor(-0.5+log10(axmin)));
        axmax = 10.^(ceil(0.5+log10(axmax)));
        loglog(resp.E_GRID,max(yeff,axmin(2)),'g-','linew',2); % epsdEG vs E
        axis([axmin(1) axmax(1) axmin(2) axmax(2)]);
        if intE,
            Gsym = 'G0';
            G_UNIT = sprintf('%s^2sr',resp.L_UNIT);
        else
            Gsym = 'G0dE';
            G_UNIT = sprintf('%s%s^2sr',resp.E_UNIT,resp.L_UNIT);
        end
        title(sprintf('%s : E0=%g (%.1f%%), %s=%g (%.1f%%)',...
            chan,results.(chan).E0,results.(chan).E0err*100,...
            Gsym,results.(chan).G0,results.(chan).G0err*100),'interpreter','none');
        xlabel(sprintf('E0, %s',resp.E_UNIT));
        ylabel(sprintf('%s, %s',Gsym,G_UNIT));
        grid on;
    end
end
