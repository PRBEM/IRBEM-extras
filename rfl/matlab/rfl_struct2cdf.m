function rfl_struct2cdf(cdfname,inst_info)
% rfl_struct2cdf(cdfname,inst_info)
% save an instrument response structure to a cdf
% Note, cell arrays of strings will be converted to rows of strings

varlist = break_down_struct('',inst_info);

extras = {};

cdfwrite(cdfname,varlist);

function varlist = break_down_struct(base,inst_info)

varlist = {};
fields = fieldnames(inst_info);
for i = 1:length(fields),
    fld = fields{i};
    if isempty(base),
        longname = fld;
    else
        longname = [base,'.',fld];
    end
    if isa(inst_info.(fld),'logical'),
        inst_info.(fld) = uint8(inst_info.(fld)); % cdfwrite doesn't understand logicals
    end
    switch(class(inst_info.(fld))),
        case 'struct',
            newvarlist = break_down_struct(longname,inst_info.(fld));
            varlist = {varlist{:},newvarlist{:}};
        case 'function_handle',
            % skip it
        case 'cell', 
            if iscellstr(inst_info.(fld)),
                varlist{end+1} = longname;
                varlist{end+1} = char(inst_info.(fld));
            else
                error('Cannot save cell array "%s" to a cdf',longname);
            end
        otherwise
            varlist{end+1} = longname;
            varlist{end+1} = inst_info.(fld);
    end
end
