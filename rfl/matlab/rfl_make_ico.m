% make a sample structure (inst_info) containing the response function
% metadata for ICO dosimeters
clear inst_info

% top-level info
inst_info.FORMAT_VERSION = '1.0.1';
inst_info.CHANNEL_NAMES = {'ELEC1','ELEC2','ELEC3','ELEC4','ELEC5'};
inst_info.L_UNIT = 'cm';
inst_info.E_UNIT = 'MeV';
inst_info.REFERENCES = {'ICO-F2 Dosimeters. PI Bern Blake'
    'Response file prepared by Paul O''Brien, paul.obrien@aero.org'
    'GEANT4/EGSnrc simulations provided by Mark Looper'
    sprintf('Created %s',datestr(now))};
inst_info.DEAD_TIME_PER_COUNT = 0;
inst_info.DEAD_TYPE = 'BLOCKING';
inst_info.COUNTS_MAX = inf; % the actual COUNTS_MAX is something else, but I don't know what.
inst_info.SPECIES = {'ELE'};
inst_info.ELE.RESP_TYPE = '[E]';
inst_info.ELE.E_TYPE = 'TBL';
inst_info.ELE.BIDIRECTIONAL = 'FALSE';
inst_info.ELE.E_GRID = [0.1,0.12,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.8,1,1.2,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,9,10,12,15,20,25,30];
% Channel-specific info
inst_info.ELEC1.ELE.G = 0.0607;
inst_info.ELEC1.ELE.E0 = 0.95;
inst_info.ELEC1.ELE.E1 = Inf;
inst_info.ELEC1.ELE.XCAL = 0.386435;
inst_info.ELEC1.ELE.XCAL_RMSE = 0.92801;
inst_info.ELEC1.ELE.EPS = [0,0,0,0,0,0,0,0,0,0.073055,0.377304,0.556567,0.756443,0.914612,0.994849,0.944417,0.910028,0.89857,0.994849,1.04068,0.985669,1.02464,0.944417,1.06591,1.07966,0.981085,0.983377,1.15759,1.06361,0.955876,0.983377,1.0361,1.04068];
inst_info.ELEC2.ELE.G = 0.0643;
inst_info.ELEC2.ELE.E0 = 1.97;
inst_info.ELEC2.ELE.E1 = Inf;
inst_info.ELEC2.ELE.XCAL = 0.275165;
inst_info.ELEC2.ELE.XCAL_RMSE = 1.02817;
inst_info.ELEC2.ELE.EPS = [0,0,0,0,0,0,0,0,0,0,0,0,0.0328898,0.312776,0.580507,0.735359,0.806254,0.884103,1.02867,1.04257,0.986968,0.853514,1.04257,0.956391,1.08983,1.02589,1.04814,0.995303,0.986968,1.05647,0.897998,1.07871,1.02589];
inst_info.ELEC3.ELE.G = 0.4263;
inst_info.ELEC3.ELE.E0 = 3.52;
inst_info.ELEC3.ELE.E1 = Inf;
inst_info.ELEC3.ELE.XCAL = 1.89545;
inst_info.ELEC3.ELE.XCAL_RMSE = 0.971318;
inst_info.ELEC3.ELE.EPS = [0,0,0,0,0,0,0,0,0,0,0,2.55588e-005,3.19494e-005,0.000242818,0.00640263,0.0940327,0.271441,0.47298,0.659428,0.758466,0.830993,0.924607,0.935455,0.962936,1.00128,1.00895,1.06197,1.04793,0.997439,1.00256,0.977628,0.942496,0.959739];
inst_info.ELEC4.ELE.G = 0.4368;
inst_info.ELEC4.ELE.E0 = 5.45;
inst_info.ELEC4.ELE.E1 = Inf;
inst_info.ELEC4.ELE.XCAL = 1.83162;
inst_info.ELEC4.ELE.XCAL_RMSE = 0.79326;
inst_info.ELEC4.ELE.EPS = [0,0,0,0,0,0,0,0,0,0,0,0,4.76644e-005,9.53287e-005,0.000200186,0.000314584,0.00100095,0.00647286,0.0481692,0.147873,0.286368,0.42354,0.571202,0.676269,0.764545,0.867966,0.973313,0.991414,0.962824,1.01621,1.03623,0.991414,0.993321];
inst_info.ELEC5.ELE.G = 0.3616;
inst_info.ELEC5.ELE.E0 = 6.75;
inst_info.ELEC5.ELE.E1 = Inf;
inst_info.ELEC5.ELE.XCAL = 1.01426;
inst_info.ELEC5.ELE.XCAL_RMSE = 0.813318;
inst_info.ELEC5.ELE.EPS = [0,0,0,0,0,0,0,0,0,0,0,0,9.07907e-005,9.07907e-005,0.000207515,0.000259394,0.000518806,0.000518806,0.0013489,0.00249022,0.0075227,0.0314532,0.0952006,0.180288,0.297011,0.400769,0.619449,0.813225,0.960434,0.940338,0.999995,0.993513,1.06615];


outpath = fileparts(which([mfilename,'.m']));
cdfname = [outpath,filesep,'ico.cdf'];

rfl_struct2cdf(cdfname,inst_info);

inst_info = rfl_load_inst_info(inst_info);

% if 0, % this code block only works on my computers (TPO)
%     db = read_HEOICO_EG('ICO');
%     E_GRID = sprintf('%g,',db.EG{i}(:,1));
%     E_GRID = E_GRID(1:(end-1)); % remove trailing ,
%     fprintf('inst_info.E_GRID = [%s];\n',E_GRID);
%     
%     for i = 1:5,
%         E_GRID = sprintf('%g,',db.EG{i}(:,1));
%         E_GRID = E_GRID(1:(end-1)); % remove trailing ,
%         G = db.G0(i);
%         Gcal = mean(db.EG{i}(db.EG{i}(:,1)>db.E0(i)*2,2)); % average of calibrated G(E) for E>2*E0
%         XCAL = G/Gcal;
%         EPS = sprintf('%g,',db.EG{i}(:,2)/Gcal);
%         EPS = EPS(1:(end-1)); % remove trailing ,
%         E0 = db.E0(i);
%         E1 = inf;
%         
%         v = sprintf('ELEC%d',i);
%         fprintf('inst_info.%s.ELE.G = %g;\n',v,G);
%         fprintf('inst_info.%s.ELE.E0 = %g;\n',v,E0);
%         fprintf('inst_info.%s.ELE.E1 = %g;\n',v,E1);
%         fprintf('inst_info.%s.ELE.XCAL = %g;\n',v,XCAL);
%         fprintf('inst_info.%s.ELE.XCAL_RMSE = %g;\n',v,db.rel_err(i));
%         fprintf('inst_info.%s.ELE.EPS = [%s];\n',v,EPS);
%         
%     end
% end
