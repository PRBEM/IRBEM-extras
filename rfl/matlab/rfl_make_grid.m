function Xgrid = rfl_make_grid(xmin,xmax,options,var)
% Xgrid = rfl_make_grid(xmin,xmax,options,var)
% creates column vector of evenly spaced points
% from xmin to xmax based on options.int_method
% and options.(['d',var]) (e.g., options.dbeta)
% for trapz int_method, includes xmin and xmax

if isfield(options,'int_method'),
    int_method = options.int_method;
else
    int_method = 'trapz'; % default
end

if isfield(options,['d',var]),
    dX = options.(['d',var]);
else
    dX = (xmin-xmax)/100;
end

switch(int_method),
    case 'trapz',
        Xgrid = linspace(xmin,xmax,round((xmax-xmin)/dX)+1)'; % include end points
    otherwise
        error('Unknown int_method "%s"',options.int_method);
end
