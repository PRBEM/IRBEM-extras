function rfl_struct2cdf(cdfname,stru)
% rfl_struct2cdf(cdfname,stru)
% save an instrument response structure to a cdf
% Note, cell arrays of strings will be converted to rows of strings

varlist = break_down_struct('',stru);

extras = {};

cdfwrite(cdfname,varlist);

function varlist = break_down_struct(base,stru)

varlist = {};
fields = fieldnames(stru);
for i = 1:length(fields),
    fld = fields{i};
    if isempty(base),
        longname = fld;
    else
        longname = [base,'.',fld];
    end
    if isa(stru.(fld),'logical'),
        stru.(fld) = uint8(stru.(fld)); % cdfwrite doesn't understand logicals
    end
    switch(class(stru.(fld))),
        case 'struct',
            newvarlist = break_down_struct(longname,stru.(fld));
            varlist = {varlist{:},newvarlist{:}};
        case 'function_handle',
            % skip it
        case 'cell', 
            if iscellstr(stru.(fld)),
                varlist{end+1} = longname;
                varlist{end+1} = char(stru.(fld));
            else
                error('Cannot save cell array "%s" to a cdf',longname);
            end
        otherwise
            varlist{end+1} = longname;
            varlist{end+1} = stru.(fld);
    end
end
