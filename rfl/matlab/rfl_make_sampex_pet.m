% make a sample structure (inst_info) containing the response function
% metadata for PET instrument
clear inst_info


stop_detector = {'P1','P1','P1','P1','P3A','P3B','P3C','P3D/P3E','P4','P5','P6/P7','P8','P8','P8','P8'};
R = [1.596,1.596, 1.596, 1.596, 1.715, 1.715,  1.715, mean([1.715 1.715]), 1.2, 1.2, mean([1.2 1.2]), 1.2, 1.2, 1.2, 1.2];
D = [5.55, 5.55, 5.55, 5.55, 6.07, 6.52, 6.97, mean([7.41 7.86]),8.49, 9.07, mean([9.75 10.48]), 11.23, 11.23, 11.23, 11.23];
% pet.half_angle = {29.9, 29.9, 29.9, 29.9, 28.6, 26.9, 25.4, mean([24.1 22.8]),18.2,17.1, mean([16.0 14.9]),14.0,14.0,14.0 ,14.0};
% pet.G = {1.792, 1.792, 1.792, 1.792, 1.752, 1.545, 1.372, mean([1.225 1.101]),0.476, 0.420, mean([0.365 0.318]), 0.278,0.278,0.278,0.278};

Nchans = 15; % number of channels
% energy limits, MeV
P1R = 1.596; % P1 radius, cm
E0 = [19.0, 21.1, 23.2, 25.2, 27.4, 37.4, 45.8, 53.0, 65.2, 70.5, 76.1, 85.1,120 , 200, 300]; % lower energy bound, MeV
E1 = [21.0, 23.2, 25.2, 27.4, 37.4, 45.8, 53.0, 65.2, 70.5, 76.1, 85.1, 120,200, 300, 500]; % upper energy bound, MeV


% top-level info
inst_info.FORMAT_VERSION = '1.1.0';
inst_info.CHANNEL_NAMES = arrayfun(@(x)sprintf('CHAN%d',x),1:Nchans,'uniform',false); % generic names, CHAN1, CHAN2, ...
inst_info.L_UNIT = 'cm';
inst_info.E_UNIT = 'MeV';
inst_info.REFERENCES = {'PET Proton Channels. PI R. Mewaldt'
    'Response file prepared by Paul O''Brien, paul.obrien@aero.org'
    'Geometry provided by Mark Looper'
    'Cook et al., PET: A Proton/Electron Telescope..., IEEE Trans. Geoscience and Remote Sensing, 31(1) May, 1993'
    sprintf('Created %s',datestr(now))};
inst_info.DEAD_TIME_PER_COUNT = 0; % PET has a livetime counter
inst_info.DEAD_TYPE = 'BLOCKING'; % unknown
inst_info.COUNTS_MAX = inf; % unknown
inst_info.SPECIES = {'PROT'};
inst_info.PROT.RESP_TYPE = '[E],[TH]'; % E independent of theta (nominally)
inst_info.PROT.E_TYPE = 'WIDE';
inst_info.PROT.BIDIRECTIONAL = 'FALSE';
inst_info.PROT.TH_TYPE = 'CYL_TELE'; % Cylindrical telescope
inst_info.PROT.CROSSCALIB = 1; % unknown
inst_info.PROT.CROSSCALIB_RMSE = log(2)/2; % unknown: assume 2 standard deviations is a factor of 2

% Channel-specific info

for i = 1:Nchans,
    v = inst_info.CHANNEL_NAMES{i};
    inst_info.(v).PROT.E0 = E0(i);
    inst_info.(v).PROT.E1 = E1(i);
    inst_info.(v).PROT.EPS = 1;
    inst_info.(v).PROT.R1 = P1R;
    inst_info.(v).PROT.R2 = R(i);
    inst_info.(v).PROT.D = D(i);
    inst_info.(v).STOP_DETECTOR = stop_detector{i};
end

outpath = fileparts(which([mfilename,'.m']));
cdfname = [outpath,filesep,'sampex_pet.cdf'];

rfl_struct2cdf(cdfname,inst_info);

inst_info = rfl_load_inst_info(inst_info);
