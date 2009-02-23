function var2cdf(cdffile,var,cell2str)
% var2cdf(cdffile,var,cell2str)
% save a variable to a CDF file
% the variable will be saved in such a way
% as to be re-readable with cdf2var and retain its structure
% also, this routine is compatible with an IDL cdf2var routine.
% if cell2str, cell arrays of strings will be converted to char with strvcat

% Global Attributes
% CreatedBy = var2cdf
%    - this CDF is one we can read/write with cdf2var and var2cdf

% Variable Attributes
% size = [ # # ...]
%    - defines size of multidimensional variables
% type = <type>
%    - indicates the type of the variable:
%      float, double, etc (numeric types)
%      char (string)
%      struct (structure)
%      cell (equivalent to an IDL collection)
% fields = 'fieldname1,fieldname2,...'
%    - string of comma-delimited field names (struct only)

if nargin < 3,
    cell2str = true; % defaults to true until IDL can handle cell arrays
end

GATTRIB = struct('CreatedBy','var2cdf','CreatedAt',datestr(now),'CreatedFrom','MATLAB');
GATTRIB.STRUCTURE_NAMES = {};

varlist = {};
clear VATTRIB
VATTRIB.size = cell(0,2);
VATTRIB.type = cell(0,2);
VATTRIB.fields = cell(0,2);
[varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,'var',var,cell2str);

cdfwrite(cdffile,varlist,'VariableAttributes',VATTRIB,'GlobalAttributes',GATTRIB,'WriteMode','overwrite','Version','2.7');

% % create file & global attributes, plus first var
% %cdfwrite(cdffile,{},'GlobalAttributes',GATTRIB,'WriteMode','overwrite','Version','2.7');
% cdfwrite(cdffile,{},'WriteMode','overwrite','Version','2.7');
% for i = 1:2:length(varlist),
%     disp(i);
%     disp(varlist{i});
%     disp(varlist(i+[0 1]));
%     cdfwrite(cdffile,varlist(i+[0 1]),'WriteMode','append');
% end
% disp('writing variable atts');
% cdfwrite(cdffile,'VariableAttributes',VATTRIB,'WriteMode','append');
%

function [varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,name,var,cell2str)

% this is where the action happens

if isnumeric(var) || isa(var,'cdfepoch'),
    VATTRIB.type(size(VATTRIB.type,1)+1,:) = {name,class(var)};
    VATTRIB.size(size(VATTRIB.size,1)+1,:) = {name,size(var)};
    varlist = {varlist{:},name,var};
elseif ischar(var) || isa(var,'function_handle'),
    if isa(var,'function_handle'),
        var = func2str(var); % convert to string, then treat as string
    end
    if isempty(var),
        warning('Substituting " "(space) for empty string "%s"',name);
        var = ' ';
    end
    VATTRIB.type(size(VATTRIB.type,1)+1,:) = {name,class(var)};
    VATTRIB.size(size(VATTRIB.size,1)+1,:) = {name,size(var)};
    varlist = {varlist{:},name,var};
elseif isstruct(var),
    GATTRIB.STRUCTURE_NAMES{end+1} = name;
    flds = fieldnames(var);
    VATTRIB.type(size(VATTRIB.type,1)+1,:) = {name,'struct'};
    VATTRIB.size(size(VATTRIB.size,1)+1,:) = {name,size(var)};
    fldstr = sprintf(',%s',flds{:});
    fldstr = fldstr(2:end); % remove leading comma
    VATTRIB.fields(size(VATTRIB.fields,1)+1,:) = {name,fldstr};
    % struct has no actual variable in CDF file
    % so just put in empty numeric array
    varlist = {varlist{:},name,nan};
    if numel(var)==1,
        for j = 1:length(flds),
            [varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,[name,'.',flds{j}],var.(flds{j}),cell2str);
        end
    else
        % must treat array of structures like cell array
        for i = 1:numel(var),
            for j = 1:length(flds),
                [varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,[name,sprintf('(%d)',i),'.',flds{j}],var.(flds{j}),cell2str);
            end
        end
    end
elseif iscell(var),
    if iscellstr(var) && cell2str,
        var = strvcat(var);
        [varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,name,var,cell2str);
    else
        GATTRIB.STRUCTURE_NAMES{end+1} = name;
        VATTRIB.type(size(VATTRIB.type,1)+1,:) = {name,'cell'};
        VATTRIB.size(size(VATTRIB.size,1)+1,:) = {name,size(var)};
        % struct has no actual variable in CDF file
        % so just put in empty numeric array
        varlist = {varlist{:},name,nan};
        % save each cell as a separate variable
        sz = size(var);
        subs = cell(1,length(sz));
        for i = 1:numel(var),
            [subs{:}] = ind2sub(sz,i); % create cell array of subscripts;
            subs_str = sprintf('%d,',subs{:});
            subs_str(end) = [];
            [varlist,VATTRIB,GATTRIB] = add_to_list(varlist,VATTRIB,GATTRIB,[name,sprintf('{%s}',subs_str)],var{i},cell2str);
        end
    end
else
    error('Unable to process variable "%s" of type "%s"',name,class(var));
end
