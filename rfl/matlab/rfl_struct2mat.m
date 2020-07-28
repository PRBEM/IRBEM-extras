function rfl_struct2mat(matname,inst_info)
% rfl_struct2mat(matname,inst_info)
% save an instrument response structure to a .mat file readable by Octave
% the saved file will contain all the fields of inst_info as separate variables
% e.g., inst_info = load(matname)
% rfl_struct2mat(matname,inst_info,version_switch)
%  pass version switch to save
%  version_switch = '-v6' by default

inst_info = filter(inst_info);

if nargin < 3,
    version_switch = '-v6';
end

% decompose struct
vars = fieldnames(inst_info);

for i = 1:length(vars),
    eval([vars{i},'=inst_info.(vars{i});']);
end

save(version_switch,matname,vars{:});

function var = filter(var)

% recursively remove function pointers

if isstruct(var),
    flds = fieldnames(var);
    for i = 1:length(flds),
        fld = flds{i};
        switch(class(var.(fld))),
            case 'function_handle',
                var = rmfield(var,fld);
            case 'struct', % convert each field
                var.(fld) = filter(var.(fld));
            otherwise
                % leave s unchanged
        end
    end
elseif iscell(var),
    for i = 1:length(var),
        var{i} = filter(var{i});
    end
end % otherwise, leave it alone (does not handle cell arrays of function handles)
