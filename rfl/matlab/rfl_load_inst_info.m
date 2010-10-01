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

% define fields that can propagate down to channel/species  responses
prop_flds = {'L_UNIT','E_UNIT','DEAD_TIME_PER_COUNT','DEAD_TYPE','COUNTS_MAX','XCAL','XCAL_RMSE',...
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
            if ~isfield(inst_info.(chan).(sp),fld) && isfield(inst_info.(sp),fld),
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
            error('rfl_load_inst_info:Error5','Missing channel property: RESP_TYPE');
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
                error('rfl_load_inst_info:Error6','Unknown RESP_TYPE: %s',inst_info.(chan).(sp).RESP_TYPE);
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
