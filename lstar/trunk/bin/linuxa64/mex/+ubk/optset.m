function [ opts ] = optset( varargin )
%OPTSET Construct option structure
%   That can be passed to ubk.* routines.
%
%   opts = optset() returns option structure skeleton.
%
%   opts = optset(param1, value1, ...) constructs option structure based on
%   the parameter-value pairs. The parameters should be a string.
%   Unrecognized parameters are ignored and emits error whose message ID
%   is 'optset:UnrecognizedParameter.'
%
%   opts = optset(opts, ...) constructs option structure based on the
%   passed options. The parameter-value pairs can be followed.

%
% $Author$
% $LastChangedDate$
% $Revision$
% $Id$
%

defaultopts = struct(...
    'IONOR',[],...
    'DS',[],...
    'DXORDR',[],...
    'DYORDPHI',[],...
    'N_PHI',[],...
    'N_THETA',[],...
    'ISCARTESIANGRID',[],...
    'M_THREADS',[],...
    'N_THREADS',[],...
    'CO_SYSTEM',[]);

opts = defaultopts;
%% Empty arguments, return available options
if ~nargin, return, end

keyvalues = varargin;
%% If the first argument is struct
if isstruct(keyvalues{1})
    newopts = keyvalues{1};
    keyvalues(1) = [];
    if length(newopts)~=1
        error('optset:InvalidArgument',...
            'struct should be a scalar.')
    end
    % Concatenate
    flds = fieldnames(newopts);
    for aFldc=flds(:)'
        aFld = aFldc{1};
        aValue = newopts.(aFld);
        AFLD = upper(aFld);
        if isfield(opts,AFLD)
            opts.(AFLD) = aValue;
        else
            error('optset:UnrecognizedParameter',...
                'Unrecognized parameter "%s".',...
                aFld)
        end
    end
end

%% Construct opt struct
if ~isempty(keyvalues)
    % Check if the parameters are even number.
    if mod(length(keyvalues),2)
        error('optset:InvalidArgument',...
            'Invalid parameter-value pairs.')
    end
    keyvalues = reshape(keyvalues, 2,length(keyvalues)/2);
    
    for aKeyValue=keyvalues
        if ~ischar(aKeyValue{1})
            error('optset:InvalidArgument',...
                'A parameter should be a string.')
        end
        % Concatenate
        aValue = aKeyValue{2};
        aFld = aKeyValue{1};
        AFLD = upper(aFld);
        if isfield(opts,AFLD)
            opts.(AFLD) = aValue;
        else
            error('optset:UnrecognizedParameter',...
                'Unrecognized parameter "%s".',...
                aFld)
        end
    end
end

end

