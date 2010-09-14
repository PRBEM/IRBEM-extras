function resp = rfl_init_Eseparable(resp)
% resp = rfl_init_Eseparable(resp)
% populate methods and extra properties for
% an [E],[...] response type
switch(resp.RESP_TYPE),
    case '[E],[TH,PH]',
        switch(resp.TP_TYPE),
            %            case 'RECT_TELE',
            %            case 'SLAB',
            %            case 'TBL',
            otherwise
                error('TP_TYPE %s not defined yet',resp.TP_TYPE);
        end
    case '[E],[TH]',
        switch(resp.TH_TYPE),
            %             case 'TBL',
            %             case 'PINHOLE',
            %             case 'CYL_TELE',
            %             case 'DISK',
            otherwise
                error('TH_TYPE %s not defined yet',resp.TH_TYPE);
        end
    case '[E]',
        resp = rfl_init_omni(resp);
    otherwise
        error('RESP_TYPE %s not defined yet',resp.RESP_TYPE);
end
switch(resp.E_TYPE),
    %     case 'TBL',
    case 'INT',
        resp = rfl_init_INT(resp);
    case 'WIDE',
        resp = rfl_init_WIDE(resp);
    case 'DIFF',
        resp = rfl_init_DIFF(resp);
    otherwise
        error('E_TYPE %s not defined yet',resp.E_TYPE);
end
if ~isfield(resp,'make_hEthetaphi'),
    resp.make_hEthetaphi = @make_hEthetaphi_Eseparable;
end
if ~isfield(resp,'make_hEalphabeta'),
    resp.make_hEalphabeta = @make_hEalphabeta_Eseparable;
end
if ~isfield(resp,'make_hEtheta'),
    error('make_hEtheta_Eseparable not defined yet');
    resp.make_hEtheta = @make_hEtheta_Eseparable;
end
if ~isfield(resp,'make_hEalphabeta'),
    error('make_hEalpha_Eseparable not defined yet');
    resp.make_hEalpha = @make_hEalpha_Eseparable;
end

function [hEthetaphi,result_code] = make_hEthetaphi_Eseparable(resp,Egrid,thetagrid,phigrid,options)
% combines results of make_hE and make_thetaphi
result_code = 1;
hE = resp.make_hE(resp,Egrid,options);
hthetaphi = resp.make_thetaphi(resp,thetagrid,phigrid,options);

hEthetaphi = merge_hE_hangles(hE,hthetaphi);

function [hEalphabeta,result_code] = make_hEalphabeta_Eseparable(resp,Egrid,alphagrid,betagrid,alpha0,beta0,phib,options)
% combines results of make_hE and make_alphabeta
result_code = 1;
hE = resp.make_hE(resp,Egrid,options);
halphabeta = resp.make_alphabeta(resp,alphagrid,betagrid,alpha0,beta0,phib,options);

hEalphabeta = merge_hE_hangles(hE,halphabeta);

function hEthetaphi = merge_hE_hangles(hE,hangles)

% hE has shape: NE x 1
% hangles has shape: N1 x N2
% hEthetaphi needs shape NE x N1 x N2

NE = length(Egrid);
[N1,N2] = size(hangles);

hangles = shiftdim(hangles,-1); % introduce singleton first dimension
hEthetaphi = repmat(hE,[1 N1 N2]).*repmat(hangles,[NE,1,1]);
