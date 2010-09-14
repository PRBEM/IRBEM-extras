function str = rfl_result_code_to_string(result_code)
% str = rfl_result_code_to_string(result_code)
% convert result_code to string
% for array of codes, returns cell array of strings

if ischar(result_code) || iscell(result_code), % no effect
    str = result_code;
    return
end

if length(result_code)>1,
    str = cell(size(result_code));
    for i = 1:length(str),
        str{i} = rc2str(result_code(i));
    end
else
    str = rc2str(result_code);    
end

function str = rc2str(rc)
switch(rc),
    case 1,
        str = 'Success';
    case 0,
        str = 'Unknown Error';
    case -1,
        str = 'Incorrect/Inconsistent Argument Size';
    case -2,
        str = 'Missing channel property: G';
    case -3,
        str = 'Missing global property: CHANNEL_NAMES';
    case -4,
        str = 'Missing channel property: SPECIES';
    case -5,
        str = 'Missing property: RESP_TYPE';
    case -6,
        str = 'Unknown RESP_TYPE';
    otherwise
        error('Unknown Result Code %d',rc);
end
