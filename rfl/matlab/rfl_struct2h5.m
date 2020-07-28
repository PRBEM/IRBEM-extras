function rfl_struct2h5(h5name,inst_info)
% rfl_struct2h5(h5name,inst_info)
% save an instrument response structure to HDF5 file
% the saved file will contain all the fields of inst_info as separate variables
% e.g., inst_info = load(h5name)

h5name = fullfile(h5name);

if exist(h5name,'file')
    delete(h5name);
end

writeVar(h5name,'',inst_info)

function writeVar(h5name,prefix,var)
% recursively write structure to file
options = {}; % {'TextEncoding','UTF-8'}; % text encoding introduced in later versions of matlab, and triggers a group creation bug

if isstruct(var)
    flds = fieldnames(var);
    anyvars = false;
    for i = 1:length(flds)
        if ~isequal(flds{i},'internal'), % skip fields called "internal"
            anyvars = true;
            writeVar(h5name,[prefix,'/',flds{i}],var.(flds{i}));
        end
    end
    if anyvars && ~isempty(prefix)
        h5writeatt(h5name,prefix,'isStruct',1,options{:});
    end
elseif iscell(var)
    for i = 1:length(var)
        writeVar(h5name,sprintf('%s/%d',prefix,i-1),var{i});
    end
    h5writeatt(h5name,prefix,'isArray',1,options{:});
elseif ischar(var)
    var = uint8(native2unicode(var,'UTF-8')); % convert to UTF-8 and store as unsigned int byte
    %fprintf('Writing %s size [%s]\n',prefix,num2str(size(var)));
    h5create(h5name,prefix,size(var),'Datatype',class(var),options{:});
    h5writeatt(h5name,prefix,'isString',1,options{:});
    h5write(h5name,prefix,var);
elseif islogical(var)
    var = uint8(var); % convert to uint8
    %fprintf('Writing %s size [%s]\n',prefix,num2str(size(var)));
    h5create(h5name,prefix,size(var),'Datatype',class(var),options{:});
    h5writeatt(h5name,prefix,'isBoolean',1,options{:});
    h5write(h5name,prefix,var);
elseif ismember(class(var),{'double', 'single','uint64', 'int64', 'uint32', 'int32', 'uint16','int16', 'uint8','int8'})
    %fprintf('Writing %s size [%s]\n',prefix,num2str(size(var)));
    h5create(h5name,prefix,size(var),'Datatype',class(var),options{:});
    h5write(h5name,prefix,var);
else
    %fprintf('Not writing %s, class %s\n',prefix,class(var));
end % otherwise, leave it alone (does not handle cell arrays of function handles)
