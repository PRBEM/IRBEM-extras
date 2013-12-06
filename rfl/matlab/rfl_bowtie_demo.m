% demonstrate rfl bowtie on differential channels

MeV = linspace(0.1,30,1000);

% top-level info
clear inst_info;
inst_info.FORMAT_VERSION = '1.1.0';
inst_info.L_UNIT = 'cm';
inst_info.E_UNIT = 'MeV';
inst_info.REFERENCES = {sprintf('Created %s',datestr(now))};
inst_info.DEAD_TIME_PER_COUNT = 0;
inst_info.DEAD_TYPE = 'BLOCKING';
inst_info.COUNTS_MAX = inf; % the actual COUNTS_MAX is something else, but I don't know what.
inst_info.SPECIES = {'ELE'};
inst_info.ELE.RESP_TYPE = '[E]'; % omni with efficiency table
inst_info.ELE.E_TYPE = 'TBL';
inst_info.ELE.BIDIRECTIONAL = 'FALSE';
inst_info.ELE.E_GRID = MeV;
% Channel-specific info
inst_info.CHANNEL_NAMES = {};
for iE = [1,3,9],
    for dE = [0.1,0.5,1],
        chan = sprintf('E%dMEV_d%dkeV',iE,dE*1e3);
        inst_info.CHANNEL_NAMES{end+1} = chan;
        inst_info.(chan).ELE.G = 1;
        inst_info.(chan).ELE.E0 = iE;
        inst_info.(chan).ELE.DE = dE;
        inst_info.(chan).ELE.CROSSCALIB = 1;
        inst_info.(chan).ELE.CROSSCALIB_RMSE = log(2)/2;
        inst_info.(chan).ELE.EPS = double(abs(MeV-iE) < dE/2);
    end
end

inst_info = rfl_load_inst_info(inst_info); % populate
results = rfl_bowtie(inst_info,'diff','ELE',2:5,[],'plot');
for i = 1:length(inst_info.CHANNEL_NAMES),
    chan = inst_info.CHANNEL_NAMES{i};
    ideal = inst_info.(chan).ELE;
    bowtie = results.(chan);
    fprintf('%20s, ideal(bowtie): E0=%2.f(%.2f), dE = %.2f(%.2f)\n',chan,ideal.E0,bowtie.E0,ideal.DE,bowtie.G0/ideal.G);
end