function inst_info = rfl_init_Eseparable_csym(inst_info)
% inst_info = rfl_init_Eseparable_csym(inst_info)
% populate methods and extra properties for
% an [E],[TH] response type
% inst_info is assumed to be a fully-initialized Eseparable sensor
% this function just overloads certain methods

inst_info.internal.hAthetaphi = @make_hAthetaphi_csym;
inst_info.internal.hAtheta = @make_hAtheta_csym;

switch(inst_info.RESP_TYPE),
    case '[E],[TH]',
        switch(inst_info.TH_TYPE),
            case 'TBL',
                inst_info.internal.A = @(inst_info,theta,phi)interpn(inst_info.TH_GRID,inst_info.A,theta,'linear',0); % zero outside grid, ignore phi
                %             case 'PINHOLE', %%
                %             case 'CYL_TELE', %%
            case 'SLAB',
                inst_info = init_single_element(inst_info,inst_info.H1*inst_info.W1);
            case 'DISK',
                inst_info = init_single_element(inst_info,pi*inst_info.R1^2);
            otherwise
                error('TH_TYPE %s not defined yet',inst_info.TH_TYPE);
        end
    case '[E]', 
        inst_info = rfl_init_omni(inst_info);
    otherwise
        error('RESP_TYPE %s not defined yet',inst_info.RESP_TYPE);
end


function inst_info = init_single_element(inst_info,area)
% initialize angle methods for single element telescopes
% (DISK or SLAB)
% for a single element with specified area

% Overload .A method (effective area)

inst_info.internal.area = area;

if strcmpi(inst_info.BIDIRECTIONAL,'TRUE'),
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.internal.area*cosd(theta);
    inst_info.internal.hA0 = area*2*pi;
else % this is the case Sullivan treats
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.internal.area*cosd(theta).*(theta<=90);
    inst_info.internal.hA0 = area*pi;
end


function inst_info = rfl_init_omni(inst_info)
% populate methods and variables for
% make_hA* part of response function with
% idealized omnidirectional channel

if strcmpi(inst_info.BIDIRECTIONAL,'TRUE'),
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.G/4/pi;
    inst_info.internal.hAthetaphi = @fullomni_hAthetaphi;
    inst_info.internal.hAtheta = @fullomni_hAtheta;
    inst_info.internal.hAalphabeta = @(inst_info,alphagrid,betagrid,alpha0,beta0,phib,options)inst_info.internal.hAthetaphi(inst_info,alphagrid,betagrid,options);
    inst_info.internal.hAalpha = @(inst_info,alphagrid,alpha0,beta0,phib,options)inst_info.internal.hAtheta(inst_info,alphagrid,options);
else % this is the case Sullivan treats
    inst_info.internal.A = @(inst_info,theta,phi)inst_info.G*(theta<=90)/2/pi;
    inst_info.internal.hAthetaphi = @halfomni_hAthetaphi;
    inst_info.internal.hAtheta = @halfomni_hAtheta;
    % keep standard csym separable hAalphabeta and hAalpha
end
inst_info.internal.hA0 = inst_info.G;

function [hAthetaphi,result_code] = fullomni_hAthetaphi(inst_info,thetagrid,phigrid,options)
result_code = 1;
hAthetaphi = inst_info.G/4/pi*rfl_make_deltas(thetagrid(:),options)*rfl_make_deltas(phigrid(:),options)';

function [hAtheta,result_code] = fullomni_hAtheta(inst_info,thetagrid,options)
result_code = 1;
hAtheta = inst_info.G/2*rfl_make_deltas(thetagrid(:),options);

function [hAthetaphi,result_code] = halfomni_hAthetaphi(inst_info,thetagrid,phigrid,options)
result_code = 1;
hAthetaphi = inst_info.G/2/pi*((thetagrid(:)<=90).*rfl_make_deltas(thetagrid(:),options))*rfl_make_deltas(phigrid(:),options)';

function [hAtheta,result_code] = halfomni_hAtheta(inst_info,thetagrid,options)
result_code = 1;
hAtheta = inst_info.G*(thetagrid(:)<=90).*rfl_make_deltas(thetagrid(:),options);

function [hAthetaphi,result_code] = make_hAthetaphi_csym(inst_info,thetagrid,phigrid,options)
dcost = rfl_make_deltas(-cosd(thetagrid),options);
dp = rfl_make_deltas(phigrid,options)*pi/180;

result_code = 1; % success
A = inst_info.internal.A(inst_info,thetagrid,0);
h = (A(:).*dcost(:))*dp(:)'; % outer product, column x row - NTHETA x NPHI matrix
hAthetaphi = h; % success, returns h

function [hAtheta,result_code] = make_hAtheta_csym(inst_info,thetagrid,options)
dcost = rfl_make_deltas(-cosd(thetagrid),options);
result_code = 1; % success
A = inst_info.internal.A(inst_info,thetagrid,0);
h = (A(:).*dcost(:))*2*pi;
hAtheta = h; % success, returns h
