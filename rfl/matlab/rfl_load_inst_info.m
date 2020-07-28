function [inst_info,result_code] = rfl_load_inst_info(FileName,FileType)
% [inst_info,result_code] = rfl_load_inst_info(FileName,FileType);
% Populates an inst_info object from file FileName.
% Supported file types are .cdf, .h5, .xml (and native, e.g., .mat)
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
    if nargin == 1,
        f = find(FileName=='.');
        ext = FileName((f(end)+1):end);
    else
        ext = replace(lower(FileType),'.','');
    end
    switch(ext),
        case 'cdf',
            inst_info = read_from_cdf(FileName);
        case 'h5',
            inst_info = read_from_h5(FileName);
        case 'mat',
            inst_info = load(FileName);
        otherwise
            error('FileType "%s" not yet supported',ext);
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
    var = info.Variables{i,1};
    f = find(var == '.');
    if isempty(f),
        last_part = var;
    else
        last_part = var((f(end)+1):end);
    end
    if ismember(last_part,{'CHANNEL_NAMES','SPECIES','REFERENCES'});
        data{i} = cellstr(data{i}); % convert rows of strings to cell arrays of strings
    end
    eval(['inst_info.',var,'=data{i};']);
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


function inst_info = read_from_h5(FileName)

info = h5info(FileName);
inst_info = struct();
for igroup = 1:length(info.Groups)
    inst_info.(h5leaf(info.Groups(igroup).Name)) = read_h5_group(FileName,info.Groups(igroup));
end

for iset = 1:length(info.Datasets)
    inst_info.(h5leaf(info.Datasets(iset).Name)) = read_h5_dataset(FileName,'',info.Datasets(iset));
end

function out = read_h5_group(FileName,info)
atts = h5atts(info);

if atts.isArray
    if isempty(info.Groups), % array of datasets
        out = cell(1,length(info.Datasets));
        for i = 1:length(out),
            out{i} = read_h5_dataset(FileName,info.Name,info.Datasets(i));
        end
    else % array of groups
        out = cell(length(info.Groups));
        for i = 1:length(out),
            out{i} = read_h5_group(FileName,info.Group(i));
        end
    end
elseif atts.isStruct
    out = struct();
    for igroup = 1:length(info.Groups)
        out.(h5leaf(info.Groups(igroup).Name)) = read_h5_group(FileName,info.Groups(igroup));
    end
    for iset = 1:length(info.Datasets)
        out.(info.Datasets(iset).Name) = read_h5_dataset(FileName,info.Name,info.Datasets(iset));
    end
end

function out = read_h5_dataset(FileName,path,dataset)
atts = h5atts(dataset);
out = h5read(FileName,[path,'/',dataset.Name]);
if atts.isBoolean,
    out = logical(out);
elseif atts.isString,
    out = char(out);
end

function atts = h5atts(info)
atts = struct('isBoolean',false,'isStruct',false,'isArray',false,'isString',false);
for iatt = 1:length(info.Attributes)
    atts.(info.Attributes(iatt).Name) = info.Attributes(iatt).Value;
end

function leaf = h5leaf(name)
% return last name after /
f = find(name == '/',1,'last');
if isempty(f),
    leaf = name;
else
    leaf = name((f+1):end);
end


