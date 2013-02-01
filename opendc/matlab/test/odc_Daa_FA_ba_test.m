% reproduce Figures 1 and 2 from
% Orlova and Sphrits, GRL, 2010

util = odc_util; % load utility functions

dipole_only = false; % only do dipole case
use_onera = true; % use onera to trace the field lines rather than loading files provided by Orlova
fig2_only = false; % only make Figure 2 (uses higher alpha0 resolution to match the published figure)
% note: if you have fig2_only = false, you'll get a good version of fig1
% but a crummy version of fig2. Set figt2_only to get a good version of
% fig2 relatively quickly at the expense of getting fig1 at all.

if fig2_only,
    R0s = 7;
else
    R0s = [4,6,7];
end
Kps = [6,2,NaN]; % NaN means dipole
Bmodel_label = {'T89, Kp=6','T89, Kp=2','Dipole'};
MLTs = [0,12]; % night, day
if fig2_only,
    alpha0_degs = 1:0.05:90; % do at very high angular resolution
else
    alpha0_degs = 1:2:90; % do at adequate angular resolution for making plots
end
MeV = 1;
ref_lats = 0:5:35;
dlat = 0.05; % 0.05 degree steps, per Orlova email
minNpts = 201;
colors = {'r','g','b','m',[0.5 0 0.5],[0.5 0 0],'k',[1 0.5 0]};

fl_path = '~/CrossTerms/Orlova/Field line points';

if use_onera,
    % read launch times
    LaunchTimes = []; % L, Kp, MLT, datenum
    fid = fopen([fl_path,'/Times.txt'],'rt');
    clear tmp
    while ~feof(fid),
        line = fgetl(fid);
        if isequal(line,-1),
            break;
        end
        if any(line=='='),
            eval(['tmp.',line]);
        else
            tmp = [];
        end
        if isfield(tmp,'matlabd'),
            LaunchTimes = [LaunchTimes; tmp.L,tmp.Kp,tmp.MLT,datenum(tmp.matlabd)];
        end
    end
    fclose(fid);
end

BBm_to_alpha = @(sign_cospa,B,Bm)acosd(sign_cospa*sqrt(1-min(B/Bm,1)));
da0da2 = @(B,Bm,Beq)max(0,(Bm./B-1)./(Bm./Beq-1));  % (da0/da)^2

fig1 = BigFig([1.5 1.5]);
fig2 = BigFig([1.5 1.5]);
for iMLT = 1:length(MLTs),
    MLT = MLTs(iMLT);
    switch(MLT),
        case 0,
            MLT_label = 'NIGHT';
            wave_model = struct('mode','R');
            wave_model.dB = @(L,MLT,maglat)0.1*(abs(maglat)<=15); % nT, confined to |maglat| < 15
            wave_model.normalization = 'Omega_e_eq';
            wave_model.omega_m = 0.35;
            wave_model.domega = 0.15;
            wave_model.omega1 = 0.05;
            wave_model.omega2 = 0.65;
            wave_model.directions = 'P'; % poleward waves
            wave_model.N0 = @(L,MLT,maglat)odc_Sheeley2001('trough',L,MLT,false); % #/cc, ignore L limits
        case 12,
            MLT_label = 'DAY';
            wave_model = struct('mode','R');
            wave_model.dB = @(L,MLT,maglat)0.1*(abs(maglat)<=35); % nT, confined to |maglat| < 35
            wave_model.normalization = 'Omega_e_eq';
            wave_model.omega_m = 0.2;
            wave_model.domega = 0.1;
            wave_model.omega1 = 0.1;
            wave_model.omega2 = 0.3;
            wave_model.directions = 'P'; % poleward waves
            wave_model.N0 = @(L,MLT,maglat)odc_Sheeley2001('trough',L,MLT,false); % #/cc, ignore L limits
        otherwise
            error('Don''t know what to do at MLT=%g',MLT);
    end
    for iR = 1:length(R0s)
        R0 = R0s(iR);
        alt_wave_model = setfield(wave_model,'dB',0.1); % for use w/ figure 2, w/o latitude limits, 100 pT
        Daa = nan(length(alpha0_degs),length(Kps));
        Daa_local = nan(length(alpha0_degs),length(Kps),length(ref_lats));
        for ialpha0_deg = 1:length(alpha0_degs),
            alpha0_deg = alpha0_degs(ialpha0_deg);
            tic;
            mirror_lat = util.dipole_mirror_latitude(alpha0_deg);
            for iBmodel = 1:length(Kps),
                Kp = Kps(iBmodel);
                toc;
                fprintf('%s,L=%d,alpha0=%g,Kp=%g\n',MLT_label,R0,alpha0_deg,Kp);
                tic;
                ba = nan;
                if isnan(Kp),
                    % dipole;
                    Npts = max(minNpts,ceil(mirror_lat/dlat)*2+1);
                    Beq = util.dipoleB(R0,0,0);
                    localf = @(XYZ,B,Bm,maglat,sign_cospa)odc_Daa_FA_local('e-',MeV,BBm_to_alpha(sign_cospa,B,Bm),R0,MLT,B,Beq,sign(maglat),wave_model).*da0da2(B,Bm,Beq);
                    alt_localf = @(XYZ,B,Bm,maglat,sign_cospa)odc_Daa_FA_local('e-',MeV,BBm_to_alpha(sign_cospa,B,Bm),R0,MLT,B,Beq,sign(maglat),alt_wave_model);
                    Bm = Beq/sind(alpha0_deg)^2;
                    if ~fig2_only,
                        ba = odc_bounce_average_dipole(R0,MLT,alpha0_deg,localf,'method','trace','Nlats',Npts);
                    end
                    for ilat = 1:length(ref_lats),
                        maglat = ref_lats(ilat);
                        Blocal = util.dipoleB(R0,maglat,0);
                        if Blocal > Bm,
                            Daa_local(ialpha0_deg,iBmodel,ilat) = 0;
                        else
                            Daa_local(ialpha0_deg,iBmodel,ilat) = (alt_localf([R0,0,0],Blocal,Bm,maglat,+1)+alt_localf([R0,0,0],Blocal,Bm,maglat,-1) ...
                                + alt_localf([R0,0,0],Blocal,Bm,-maglat,+1)+alt_localf([R0,0,0],Blocal,Bm,-maglat,-1))/4;
                        end
                    end
                elseif ~dipole_only, % T89
                    
                    if use_onera,
                        
                        itime = find((LaunchTimes(:,1)==R0) & (LaunchTimes(:,2)==Kp) & (LaunchTimes(:,3)==MLT));
                        launchtime = LaunchTimes(itime,4);
                        X0GSM = -cos(MLT*pi/12)*R0;
                        [trace.Lm,trace.Blocal,trace.Bmin,trace.J,trace.POSIT] = onera_desp_lib_trace_field_line('T89',[],'GSM',launchtime,X0GSM,0,0,onera_desp_lib_maginputs(Kp),0.95);
                        if trace.POSIT(1,3)>trace.POSIT(end,3), % reverse, south to north
                            trace.Blocal = trace.Blocal(end:-1:1,:);
                            trace.POSIT = trace.POSIT(end:-1:1,:);
                        end
                        trace.GSM = onera_desp_lib_coord_trans(trace.POSIT,'GEO2GSM',launchtime);
                        trace.GSM(:,2) = 0; % remove Y component
                        [trace.Bmin,trace.iBmin] = min(trace.Blocal);
                        trace.rBmin = trace.GSM(trace.iBmin,:);
                        trace.rBmin = trace.rBmin/sqrt(trace.rBmin*trace.rBmin');
                        trace.lat = acosd(min(1,trace.GSM*trace.rBmin'./sqrt(sum(trace.POSIT.^2,2)))); % angle between radius to Bmin and each point
                        trace.lat(1:trace.iBmin) = -trace.lat(1:trace.iBmin); % southern hemisphere
                        clear line
                        line.r = sqrt(sum(trace.GSM.^2,2));
                        line.lat = trace.lat;
                        line.B = trace.Blocal;
                        line.hemi = ones(size(line.r));
                        line.hemi(1:(trace.iBmin)) = -1;
                    else
                        % read northern and southern hemisphere points provided
                        % by Orlova
                        % Field_line_NH_L=4_Kp=2_MLT=0.mat
                        fl_file = sprintf('%s/Field_line_NH_L=%d_Kp=%d_MLT=%d.mat',fl_path,R0,Kp,MLT);
                        NH = load(fl_file);
                        fl_file = sprintf('%s/Field_line_SH_L=%d_Kp=%d_MLT=%d.mat',fl_path,R0,Kp,MLT);
                        SH = load(fl_file);
                        % merge
                        field_line = cat(1,SH.Field_line_SH(end:-1:2,:),NH.Field_line_NH);
                        % R (Re), latitude (radians), B (G)
                        clear line
                        line.r = field_line(:,1);
                        line.B = field_line(:,3)*1e5; % G to nT
                        line.lat = field_line(:,2)*180/pi; % latitude (not in the B/B0 sense)
                        line.hemi = ones(size(line.r));
                        line.hemi(1:(length(SH.Field_line_SH)-1)) = -1; % southern hemisphere
                    end
                    
                    
                    Beq = min(line.B);
                    Bm = Beq/sind(alpha0_deg)^2;
                    if Bm > min(line.B([1,end])),
                        ba = 0; % particle is in local loss cone
                    else
                        line.maglat = util.BB0toMagLat(line.B/Beq).*line.hemi;
                        line.GSM = [-cos(MLT*pi/12)*line.r.*cosd(line.lat),line.r*0,line.r.*sind(line.lat)]; % field line lies on X-GSM
                        
                        % interpolate onto imaglat
                        latlimS = max(-mirror_lat,line.maglat(1));
                        latlimN = min(mirror_lat,line.maglat(end));
                        Npts = max(floor(minNpts/2),ceil((latlimN-latlimS)/dlat/2))*2+1;
                        imaglat = linspace(latlimS,latlimN,Npts)';
                        B = interp1(line.maglat,line.B,imaglat,'linear');
                        GSM = interp1(line.maglat,line.GSM,imaglat,'linear');
                        lat = interp1(line.maglat,line.lat,imaglat,'linear'); % use angle in physical space instead of B/B0
                        % set up local function
                        localf = @(XYZ,B,Bm,maglat,sign_cospa)odc_Daa_FA_local('e-',MeV,BBm_to_alpha(sign_cospa,B,Bm),R0,MLT,B,Beq,sign(maglat),wave_model,'maglat',maglat).*da0da2(B,Bm,Beq);
                        % bounce average
                        if ~fig2_only,
                            ba = odc_bounce_average_trace(GSM,B,localf,sign(imaglat),'maglat',lat,'Bm',Bm);
                            assert(isreal(ba));
                        end
                        alt_localf = @(XYZ,B,Bm,maglat,sign_cospa)odc_Daa_FA_local('e-',MeV,BBm_to_alpha(sign_cospa,B,Bm),R0,MLT,B,Beq,sign(maglat),alt_wave_model,'maglat',maglat);
                        for ilat = 1:length(ref_lats),
                            lat = ref_lats(ilat);
                            Blocal = interp1(line.lat,line.B,lat,'linear');
                            maglat = util.BB0toMagLat(Blocal/Beq); % positive
                            if Blocal > Bm,
                                Daa_local(ialpha0_deg,iBmodel,ilat) = 0;
                            else
                                Daa_local(ialpha0_deg,iBmodel,ilat) = (alt_localf([R0,0,0],Blocal,Bm,maglat,+1)+alt_localf([R0,0,0],Blocal,Bm,maglat,-1) ...
                                    + alt_localf([R0,0,0],Blocal,Bm,-maglat,+1)+alt_localf([R0,0,0],Blocal,Bm,-maglat,-1))/4;
                            end
                        end
                    end
                end
                Daa(ialpha0_deg,iBmodel) = ba;
            end
        end
        figure(fig1);
        subplot(length(R0s),length(MLTs),sub2ind([length(MLTs),length(R0s)],iMLT,iR));
        h = semilogy(alpha0_degs,Daa);
        set(h(1),'color','r');
        set(h(2),'color','b');
        set(h(3),'color','k');
        set(gca,'xlim',[0,90],'xtick',0:20:90);
        if (iR==1) && (iMLT==1),
            legend(Bmodel_label{:},'location','nw');
        end
        if iMLT == 1,
            set(gca,'ylim',[1e-8,1e-2],'ytick',10.^[-8:2:-2]);
            ylabel('D_{\alpha0\alpha0},local,1/sec');
        else
            set(gca,'ylim',[1e-7,1e-2],'ytick',10.^[-7:2:-2]);
        end
        title(sprintf('%s,R_0=%g',MLT_label,R0));
        if iR == length(R0s),
            xlabel('\\alpha_0, deg');
        end
        drawnow;
    end % for iR
    figure(fig2);
    for iBmodel = 1:length(Kps),
        irow = length(Kps)-iBmodel+1;
        iR0 = find(R0s==7);
        R0 = R0s(iR);
        subplot(length(Kps),length(MLTs),sub2ind([length(MLTs),length(Kps)],iMLT,irow));
        for ilat = 1:length(ref_lats),
            plot(alpha0_degs,Daa_local(:,iBmodel,ilat),'-','color',colors{ilat},'linew',2);
            hold on;
        end
        if iMLT==1,
            axis([0 90 0 0.015]);
            ylabel('D_{\alpha\alpha},local,1/sec');
        else
            axis([0 90 0 0.01]);
        end
        set(gca,'xtick',0:20:90);
        if irow==1,
            title(sprintf('%s,R_0=%g,%s',MLT_label,R0,Bmodel_label{iBmodel}))
        else
            title(Bmodel_label{iBmodel});
        end
        if irow == length(Kps),
            xlabel('\alpha_0,(^o)');
        end
    end
end % for iMLT


figure(fig2);
subplot(length(Kps),length(MLTs),length(Kps)*length(MLTs));
p = get(gca,'pos');
leg = arrayfun(@(x)sprintf('\\lambda=%g^o',x),ref_lats,'uniform',false);
hl = legend(leg,'location','so','orientation','horiz');
set(gca,'pos',p);
set(hl,'pos',[0.1 0.0 0.8 0.025]);
