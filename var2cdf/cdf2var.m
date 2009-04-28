function var = cdf2var(cdffile,preserve_format)
% var = cdf2var(cdffile)
% var = cdf2var(cdffile,preserve_format)
% load a variable from a CDF file created by var2cdf
% only files saved by the IDL or Matlab version of var2cdf
% will work because they include special metadata that
% allows the routine to reconstruct structured variables.
%
% preserve_format (default = false)
%   if true, preserves the numeric format in the CDF.
%   if false, converts all numeric format data into double precision

if nargin < 2,
    preserve_format = false;
end


info = cdfinfo(cdffile);
rawvars = info.Variables(:,1); % list of variable names
[vars,isort] = sort(rawvars); % sort list to ensure struct appears before sub-structs

var = [];
data = cdfread(cdffile,'ConvertEpochToDatenum', true); % all variables, in rawvars order
for ivar = 1:length(vars),
    if iscell(data{isort(ivar)}),
        data{isort(ivar)} = cell2mat(data{isort(ivar)});
    end
    if isfield(info.VariableAttributes,'size'),
        isize = strmatch(vars{ivar},info.VariableAttributes.size(:,1),'exact'); % find var's index in data cell array
        sz = info.VariableAttributes.size{isize,2};
        if numel(sz)==1,
            sz = [sz 1];
        end
    else
        sz = []; % default is empty
    end
    if isfield(info.VariableAttributes,'type'),
        itype = strmatch(vars{ivar},info.VariableAttributes.type(:,1),'exact'); % find var's index in data cell array
        type_str = lower(info.VariableAttributes.type{itype,2});
    else
        if isfield(info.GlobalAttributes,'STRUCTURE_NAMES') && ismember(vars{ivar},info.GlobalAttributes.STRUCTURE_NAMES),
            type_str = 'struct';
        else
            type_str = class(data{isort(ivar)});
        end
    end
    switch(type_str),
        case {'cell'},
            % set up empty cell
            cmd = sprintf('%s=cell([%s]);',vars{ivar},sprintf('%d ',sz));
            eval(cmd); % execute
        case {'struct'},
            % do nothing, for now
        case {'char','uchar',... % character types (some of these strings are guesses)
                'byte','int1','int2','int4','uint1','unit2','uint4',... % CDF int types
                'int8','int16','int32','int64','uint8','uint16','uint32','uint64', ...
                'logical',...
                'real4','real8','float','double','single',... % real types
                'unknown'}, % unknown, for incoming idl variables
            vardata = data{isort(ivar)};
            if iscell(vardata),
                vardata = cell2mat(vardata);
            end
            if isnumeric(vardata),
                if ~preserve_format,
                    vardata = double(vardata); % force numeric data to Matlab's standard format
                end
                if ~isempty(sz) && ~isempty(vardata),
                    vardata = reshape(vardata,sz);
                end
            end
            cmd = [vars{ivar},'=','vardata;'];
            eval(cmd); % copy data;
        case {'cdf_epoch', 'cdf_epoch16'} %  or is it 'epoch','epoch16'?
            error('Don''t know how to handle type "%s"',info.VariableAttributes.type{itype,2});
        otherwise
            error('Don''t know how to handle type "%s"',info.VariableAttributes.type{itype,2});
    end
end
return;
