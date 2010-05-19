% Example script to find conjunctions between a HEO and GPS vehicle
% This examples uses made-up ephemeris, but in a real application
% one would load the ephemeris and magnetic coordinates from data files.
% Because of the long times it takes to compute magnetic coordinates,
% this example is limited to one month. The kdtree algorithm is actually
% fast enough that one can operate on years of data at a time.

%% create fake ephemeris for one month
start_date = datenum(2000,1,1);
end_date = datenum(2000,2,1);
sysaxes = 'GEI';
% HEO
fprintf('Setting up simulated HEO ephemeris and magnetic coordinates\n....(this may take a while)\n');
HEO.deltasec = 60; % simulate 1-minute data
HEO.elements = struct('type','solar','i',asin(sqrt(4/5))*180/pi,'A_p',2195,'A_a',38268,'H_a',0,'H_i',0,'M0',180);
HEO.pos = onera_desp_lib_sgp4_ele(HEO.elements,start_date,end_date,HEO.deltasec,sysaxes);
[HEO.Lm,HEO.Lstar,HEO.Blocal,HEO.Bmin,HEO.J,HEO.MLT] = onera_desp_lib_make_lstar('opq','',sysaxes,HEO.pos.date,HEO.pos.X(:,1),HEO.pos.X(:,2),HEO.pos.X(:,3));

% GPS
fprintf('Setting up simulated GPS ephemeris and magnetic coordinates\n....(this may take a while)\n');
GPS.deltasec = 5*60; % simulate 5 minute data
GPS.elements = struct('type','solar','i',55,'A_p',20200,'A_a',20200,'H_a',0,'H_i',12,'M0',270);
GPS.pos = onera_desp_lib_sgp4_ele(GPS.elements,start_date,end_date,GPS.deltasec,sysaxes);
[GPS.Lm,GPS.Lstar,GPS.Blocal,GPS.Bmin,GPS.J,GPS.MLT] = onera_desp_lib_make_lstar('opq','',sysaxes,GPS.pos.date,GPS.pos.X(:,1),GPS.pos.X(:,2),GPS.pos.X(:,3));


%% Set up Sebastien's requirements
fprintf('Setting up conjunction requirements from Friedel et al., 2005\n');
% in Friedel et al., SPACE WEATHER, VOL. 3, S09B04,
% doi:10.1029/2005SW000153, 2005
% L < 6
maxL = 6;
% within 2 hours of 0600 or 1800
dMLT = 2;
MLT1 = 6;
MLT2 = 18;

% delta L < 0.1
deltaL = 0.1;
% deltaBB0 < 0.1
deltaBB0 = 0.1;
% deltaT = 2 hours
deltaT = 2/24; % time unit is days
% my requirement: deltaMLT < 4, raise limit on Kp
deltaMLT = 4;

% ignore Kp example. Kp restrictions should be handled in pre-processing,
% i.e., by selecting only points from each spacecraft time series that
% satisfy the Kp<2 over preceding 48 hours condition, or some more relaxed
% version thereof.

%% define coordinates and distance scale
fprintf('Defining coordinates, distance scale, and 4-d search space\n');
%X: time, L, MLT, B/B0
coords = {'time','L','MLT','B/B0'};
% set distance scale to maximum allowed delta for each coordinate
DistScale = [1./deltaT, 1./deltaL, 1./deltaMLT, 1./deltaBB0];

%% set up databases of points in 4-D space for conjunctions
HEO.X = [HEO.pos.date,abs(HEO.Lm),abs(HEO.MLT),HEO.Blocal./HEO.Bmin];
GPS.X = [GPS.pos.date,abs(GPS.Lm),abs(GPS.MLT),GPS.Blocal./GPS.Bmin];

%% filter for points that are in limits and remember their indices in
% original data set
fprintf('Filtering 4-D coordinates for valid points, valid L, and valid MLT\n');
inlimits = @(db)all(isfinite(db.X),2) & (db.X(:,2)<maxL) & ...
    ((abs(db.X(:,3)-MLT1)<dMLT) | (abs(db.X(:,3)-MLT2)<dMLT));
HEO.finlimits = find(inlimits(HEO));
HEO.X = HEO.X(HEO.finlimits,:);
GPS.finlimits = find(inlimits(GPS));
GPS.X = GPS.X(GPS.finlimits,:);

%% set up kdtree for fast search
fprintf('Setting up Kdtree\n');
tic;
tree = kdtree('build',HEO.X);
toc;

%% find nearest neighbors
fprintf('Finding Nearest Neighbors\n');
tic;
% find heo points nearest each gps point
k = 1; % number of nearest neighbors - only need first, i.e., best conjunction
% use "max" distance metric, which sets R2 to the square of the largest
% difference in X between the nearest neighbor
[index,R2,XNN] = kdtree('multikNN',HEO.X,tree,GPS.X,k,DistScale,'max');
iconj = find(R2(1,:)'<1); % nearest point meets all criteria because it's max difference is < 1
iXconj = index(1,iconj); % get index of conjunction point into HEO.X
iHEOconj = HEO.finlimits(iXconj); % convert into index into HEO.Lm, HEO.pos, etc.
iGPSconj = GPS.finlimits(iconj); % convert index into GPS.Lm, GPS.pos, etc.
fprintf('%d conjunctions found\n',length(iconj));
toc;


%% plot conjuctions, one subplot per coordinate
fprintf('Plotting\n');
figure;
for i = 1:4,
    subplot(2,2,i);
    scatter(HEO.X(iXconj,i),GPS.X(iconj,i),4,1:length(iconj),'linew',3);
    ax = axis;
    if strcmpi(coords{i},'MLT'),
        ax = [0 24 0 24];
    else
        ax([1,3]) = min(ax([1,3]));
        ax([2,4]) = max(ax([2,4]));
    end
    axis(ax);
    grid on;
    if strcmpi(coords{i},'time'),
        datetick('x');
        datetick('y');
        title(sprintf('%d conjuntions',length(iconj)));
    else
        title(coords{i});
    end
    if i>= 3,
        xlabel('HEO');
    end
    if rem(i,2),
        ylabel('GPS');
    end
end
