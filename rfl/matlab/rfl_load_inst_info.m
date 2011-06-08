function [inst_info,result_code] = rfl_load_inst_info(FileName,FileType)
% [inst_info,result_code] = rfl_load_inst_info(FileName,FileType);
% Populates an inst_info object from file FileName.
% Supported file types are .cdf, .hd5, .xml (and native, e.g., .mat)
% FileType can be specified with or without the .
% if FileType is omitted, the extension of FileName will be used to
% determine what file type it is.
% FileName can also be a structure, in which case
% the returned inst_info will be "fleshed out" with
% appropriate function pointers (methods),
% and all properties propagated to the channels from their parent

result_code = 1; % assume success! (error handling isn't really implemented yet)

if isstruct(FileName),
    inst_info = FileName;
else
    f = find(FileName=='.');
    ext = FileName((f(end)+1):end);
    switch(lower(ext)),
        case 'cdf',
            inst_info = read_from_cdf(FileName);
        otherwise
            error('Extension "%s" not yet supported',ext);
    end
end

if ~isfield(inst_info,'CHANNEL_NAMES'),
    error('Missing global property: CHANNEL_NAMES');
end

renames = struct('XCAL','CROSSCALIB','XCAL_RMSE','CROSSCALIB_RMSE');
% XCAL to CROSSCALIB rename following upgrade from RFL v1.0.0 to v1.1.0
%  in response to change from format specification v1.0.1 to v1.1.0

inst_info = recursive_rename(inst_info,renames);


% define fields that can propagate down to channel/species  responses
prop_flds = {'L_UNIT','E_UNIT','DEAD_TIME_PER_COUNT','DEAD_TYPE','COUNTS_MAX','CROSSCALIB','CROSSCALIB_RMSE',...
    'RESP_TYPE','ETP_TYPE','ET_TYPE','E_TYPE','TP_TYPE','TH_TYPE','E_GRID','TH_GRID','PH_GRID','EPS',...
    'R','A','G','E0','E1','DE','R1','R2','W1','W2','H1','H2','D','BIDIRECTIONAL'};


for ichan = 1:length(inst_info.CHANNEL_NAMES),
    chan = inst_info.CHANNEL_NAMES{ichan};
    % copy inst_info to chan
    if ~isfield(inst_info.(chan),'SPECIES'),
        inst_info.(chan).SPECIES = inst_info.SPECIES;
    end
    if ~isfield(inst_info,'SPECIES'),
        error('Missing channel property: SPECIES');
    end
    for ifld = 1:length(prop_flds),
        fld = prop_flds{ifld};
        if ~isfield(inst_info.(chan),fld) && isfield(inst_info,fld),
            inst_info.(chan).(fld) = inst_info.(fld);
        end
    end
    % copy chan to species
    for isp = 1:length(inst_info.(chan).SPECIES),
        sp = inst_info.(chan).SPECIES{isp};
        for ifld = 1:length(prop_flds),
            fld = prop_flds{ifld};
            if ~isfield(inst_info.(chan).(sp),fld) && isfield(inst_info.(chan),fld),
                inst_info.(chan).(sp).(fld) = inst_info.(chan).(fld);
            end
            % copy inst_info.(sp) to inst_info.(chan).(sp)
            if ~isfield(inst_info.(chan).(sp),fld) && isfield(inst_info,sp) && isfield(inst_info.(sp),fld),
                inst_info.(chan).(sp).(fld) = inst_info.(sp).(fld);
            end
        end
    end
end

for ichan = 1:length(inst_info.CHANNEL_NAMES),
    chan = inst_info.CHANNEL_NAMES{ichan};
    for isp = 1:length(inst_info.(chan).SPECIES),
        sp = inst_info.(chan).SPECIES{isp};
        if ~isfield(inst_info.(chan).(sp),'RESP_TYPE'),
            error('Missing channel property: RESP_TYPE');
        end
        % initialize to inseparable, then let other initializations overload methods
        inst_info.(chan).(sp) = rfl_init_inseparable(inst_info.(chan).(sp));
        switch(inst_info.(chan).(sp).RESP_TYPE),
            case {'[E,TH,PH]'},
                % do nothing, already fully initialized
            case {'[E,TH]'},
                inst_info.(chan).(sp) = rfl_init_inseparable_csym(inst_info.(chan).(sp));
            case {'[E]','[E],[TH]','[E],[TH,PH]'}
                inst_info.(chan).(sp) = rfl_init_Eseparable(inst_info.(chan).(sp));
            otherwise
                error('Unknown RESP_TYPE: %s',inst_info.(chan).(sp).RESP_TYPE);
        end
    end
end

function inst_info = read_from_cdf(FileName)
% inst_info = read_from_cdf(FileName)
% load inst_info from a CDF
inst_info = [];

info = cdfinfo(FileName);
data = cdfread(FileName);
for i = 1:length(info.Variables),
    if ismember(info.Variables{i,1},{'CHANNEL_NAMES','SPECIES','REFERENCES'}),
        data{i} = cellstr(data{i}); % convert rows of strings to cell arrays of strings
    end
    eval(['inst_info.',info.Variables{i,1},'=data{i};']);
end

function var = recursive_rename(var,renames)
% recursively rename fields of struct var
% following mapping in struct renames

flds = fieldnames(renames);
if isstruct(var),
    for i = 1:length(flds),
        fld = flds{i};
        if isfield(var,fld),
            var.(renames.(fld)) = var.(fld);
            var = rmfield(var,fld);
        end
    end
elseif iscell(var),
    for i = 1:length(var),
        var{i} = recursive_rename(var{i});
    end
end
