% make one instance of each type of response function

clear inst_info

%% top-level info
inst_info.FORMAT_VERSION = '1.1.0';
inst_info.L_UNIT = 'cm';
inst_info.E_UNIT = 'MeV';
inst_info.REFERENCES = {'Made-up sensors of all types'
    'Response file prepared by Paul O''Brien, paul.obrien@aero.org'
    sprintf('Created %s',datestr(now))};
inst_info.DEAD_TIME_PER_COUNT = 0; % unknown
inst_info.DEAD_TYPE = 'BLOCKING'; % unknown
inst_info.COUNTS_MAX = inf; % unknown
inst_info.SPECIES = {'PROT'};
inst_info.CROSSCALIB = 1; % unknown
inst_info.CROSSCALIB_RMSE = log(2)/2; % unknown: assume 2 standard deviations is a factor of 2

clear tmp_chans


%% [E,TH,PH] - arbitrary geometry
% "ARB" ETP_TYPE TBL
tmp_chans.ARB.RESP_TYPE = '[E,TH,PH]'; % E response depends on theta, phi
tmp_chans.ARB.ETP_TYPE = 'TBL';
tmp_chans.ARB.E_GRID = [1 2 3 4 5];
tmp_chans.ARB.TH_GRID = (0:1:15)'; % response out to 15 degrees
tmp_chans.ARB.PH_GRID = (0:1:360)'; % response at all angles
[EI,THI,PHI] = ndgrid(tmp_chans.ARB.E_GRID,tmp_chans.ARB.TH_GRID,tmp_chans.ARB.PH_GRID);
tmp_chans.ARB.R = 3*max(0,atan((EI-2)*10)).*exp(-3*(THI/15).^2).*(1+cosd(PHI).^2);


%% [E,TH] - arbitrary geometry, cylindrically symmetric
% "ARBC" ET_TYPE TBL
tmp_chans.ARBC.RESP_TYPE = '[E,TH]'; % E response depends on theta, cylindrically symmetric
tmp_chans.ARBC.ET_TYPE = 'TBL';
tmp_chans.ARBC.E_GRID = [1 2 3 4 5];
tmp_chans.ARBC.TH_GRID = (0:1:15)'; % response out to 15 degrees
[EI,THI] = ndgrid(tmp_chans.ARBC.E_GRID,tmp_chans.ARBC.TH_GRID);
tmp_chans.ARBC.R = 3*max(0,atan((EI-2)*10)).*exp(-3*(THI/15).^2);

%% [E],[TH,PH] - separable energy/angle response
% "ESEP" E_TYPE INT, TP_TYPE TBL
tmp_chans.ESEP.RESP_TYPE = '[E],[TH,PH]'; % E independent of theta
tmp_chans.ESEP.E_TYPE = 'INT';
tmp_chans.ESEP.E0 = 20;
tmp_chans.ESEP.EPS = 1;
tmp_chans.ESEP.TP_TYPE = 'TBL'; % table of theta,phi response
tmp_chans.ESEP.TH_GRID = (0:1:15)'; % response out to 15 degrees
tmp_chans.ESEP.PH_GRID = (0:1:360)'; % response at all angles
[THI,PHI] = ndgrid(tmp_chans.ESEP.TH_GRID,tmp_chans.ESEP.PH_GRID);
tmp_chans.ESEP.A = exp(-3*(THI/15).^2).*(1+cosd(PHI).^2);

% "RECT_TELE" E_TYPE WIDE, TP_TYPE RECT_TELE, BIDIRECTIONAL FALSE
tmp_chans.RECT_TELE.RESP_TYPE = '[E],[TH,PH]'; % E independent of theta
tmp_chans.RECT_TELE.E_TYPE = 'WIDE';
tmp_chans.RECT_TELE.E0 = 20;
tmp_chans.RECT_TELE.E1 = 22;
tmp_chans.RECT_TELE.EPS = 1;
tmp_chans.RECT_TELE.TP_TYPE = 'RECT_TELE'; % two-element rectangular telescope
tmp_chans.RECT_TELE.W1 = 1;
tmp_chans.RECT_TELE.H1 = 2;
tmp_chans.RECT_TELE.W2 = 3;
tmp_chans.RECT_TELE.H2 = 4;
tmp_chans.RECT_TELE.D = 5;
tmp_chans.RECT_TELE.BIDIRECTIONAL = 'FALSE';

%% [E],[TH] - separable energy/angle response, cylindrically symmetric
% "PINHOLE" E_TYPE DIFF, TH_TYPE PINHOLE, BIDIRECTIONAL FALSE
tmp_chans.PINHOLE.RESP_TYPE = '[E],[TH]'; % E independent of theta, cylindrically symmetric
tmp_chans.PINHOLE.E_TYPE = 'WIDE';
tmp_chans.PINHOLE.E0 = 20;
tmp_chans.PINHOLE.E1 = 22;
tmp_chans.PINHOLE.EPS = 1;
tmp_chans.PINHOLE.TH_TYPE = 'PINHOLE'; % pinhole angular response
tmp_chans.PINHOLE.G = 1;
tmp_chans.PINHOLE.BIDIRECTIONAL = 'FALSE';

% "CYL_TELE" E_TYPE WIDE, TH_TYPE CYL_TELE, BIDIRECTIONAL FALSE
tmp_chans.CYL_TELE.RESP_TYPE = '[E],[TH]'; % E independent of theta, cylindrically symmetric
tmp_chans.CYL_TELE.E_TYPE = 'WIDE';
tmp_chans.CYL_TELE.E0 = 20;
tmp_chans.CYL_TELE.E1 = 22;
tmp_chans.CYL_TELE.EPS = 1;
tmp_chans.CYL_TELE.TH_TYPE = 'CYL_TELE'; % Cylindrical telescope
tmp_chans.CYL_TELE.R1 = 1;
tmp_chans.CYL_TELE.R2 = 2;
tmp_chans.CYL_TELE.D = 6;
tmp_chans.CYL_TELE.BIDIRECTIONAL = 'FALSE';

% "DISK" E_TYPE WIDE, TH_TYPE DISK, BIDIRECTIONAL FALSE
tmp_chans.DISK.RESP_TYPE = '[E],[TH]'; % E independent of theta, cylindrically symmetric
tmp_chans.DISK.E_TYPE = 'WIDE';
tmp_chans.DISK.E0 = 20;
tmp_chans.DISK.E1 = 22;
tmp_chans.DISK.EPS = 1;
tmp_chans.DISK.TH_TYPE = 'DISK'; % single-element disk
tmp_chans.DISK.R1 = 1;
tmp_chans.DISK.BIDIRECTIONAL = 'FALSE';

% "SLAB" E_TYPE WIDE, TH_TYPE SLAB, BIDIRECTIONAL FALSE
tmp_chans.SLAB.RESP_TYPE = '[E],[TH]'; % E independent of theta, cylindrically symmetric
tmp_chans.SLAB.E_TYPE = 'WIDE';
tmp_chans.SLAB.E0 = 20;
tmp_chans.SLAB.E1 = 22;
tmp_chans.SLAB.EPS = 1;
tmp_chans.SLAB.TH_TYPE = 'SLAB'; % single-element slab
tmp_chans.SLAB.W1 = 1;
tmp_chans.SLAB.H1 = 2;
tmp_chans.SLAB.BIDIRECTIONAL = 'FALSE';

%% [E] - omnidirectional response
% "OMNI" E_TYPE TBL, BIDIRECTIONAL FALSE
tmp_chans.OMNI.RESP_TYPE = '[E]'; % E independent of angle
tmp_chans.OMNI.E_TYPE = 'TBL';
tmp_chans.OMNI.E_GRID = [1 2 3 4 5];
tmp_chans.OMNI.EPS =    [0 0.5 1 1 1 ];
tmp_chans.OMNI.G = 3;
tmp_chans.OMNI.BIDIRECTIONAL = 'FALSE';

%% make bidirectional versions of all where appropriate

names = fieldnames(tmp_chans);
for i = 1:length(names),
    v = names{i};
    if isfield(tmp_chans.(v),'BIDIRECTIONAL'),
        vb = [v,'B'];
        resp = tmp_chans.(v);
        resp.BIDIRECTIONAL = 'TRUE';
        tmp_chans.(vb) = resp;
    end
end

% set channel names
inst_info.CHANNEL_NAMES = sort(fieldnames(tmp_chans));
for i = 1:length(inst_info.CHANNEL_NAMES),
    v = inst_info.CHANNEL_NAMES{i};
    inst_info.(v).PROT = tmp_chans.(v);
end

outpath = fileparts(which([mfilename,'.m']));
cdfname = [outpath,filesep,'all_types.cdf'];

rfl_struct2cdf(cdfname,inst_info);
rfl_struct2mat(strrep(cdfname,'.cdf','.mat'),inst_info);

inst_info = rfl_load_inst_info(inst_info);

for i = 1:length(inst_info.CHANNEL_NAMES),
    v = inst_info.CHANNEL_NAMES{i};
    sp = 'PROT';
    resp = inst_info.(v).(sp);
    if isfield(resp,'internal') && isfield(resp.internal,'hA0'),
        G = resp.internal.hA0;
    elseif isfield(resp,'G'),
        G = resp.G;
    elseif isfield(resp,'R'), % guess geometric factor from RE(E)
        switch(resp.RESP_TYPE),
            case '[E,TH]',
                G = 2*pi*max(trapz(-cosd(resp.TH_GRID),resp.R,2));
            case '[E,TH,PH]',
                G = max(trapz(resp.PH_GRID*pi/180,trapz(-cosd(resp.TH_GRID),resp.R,2),3));
            otherwise
                error('Huh? Cannot estimate G for resp type %s',resp.RESP_TYPE);
        end
    end
    fprintf('%s.%s.Geff = %g %s^2 sr\n',v,sp,G,resp.L_UNIT);
end
