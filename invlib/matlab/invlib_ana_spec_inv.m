% demo ana_spec_inv - just like specinv_test.c

infile = 'specinv_test.in1';

if ~exist(infile,'file'),
    error('Unable to locate "%s". Please copy it into matlab path or edit "%s.m"',infile,mfilename);
end

fid = fopen(infile,'r');
NC = fscanf(fid,'%li',1);
NE = fscanf(fid,'%li',1);
NEout = fscanf(fid,'%li',1);

c = nan(NC,1);
dc = nan(NC,1);
Egrid = nan(NE,1);
H = nan(NC*NE,1);
b = nan(NC,1);
Eout = nan(NEout,1);
flux = nan(NEout,1);
dlogflux = nan(NEout,1);

c(:) = fscanf(fid,'%lf',NC);
dc(:) = fscanf(fid,'%lf',NC);
Egrid(:) = fscanf(fid,'%lf',NE);
H(:) = fscanf(fid,'%lf',NC*NE);
b(:) = fscanf(fid,'%lf',NC);
Eout(:) = fscanf(fid,'%lf',NEout);
fclose(fid);

int_params = int32(zeros(7,1));
real_params = nan(1,1);

int_params(1+0) = NC;
int_params(1+1) = NE;
int_params(1+2) = NEout;
int_params(1+3) = 1+2; % analtyical functions bitmap*/
int_params(1+4) = 0; % minimizer, 0=BFGS, 3=NM */
int_params(1+5) = 1000; % maximumn # of iterations */
int_params(1+6) = 0; % 0 = verbose off, no text output */

real_params(1+0) = 0.511; % electron rest energy, MeV */

if ~libisloaded('invlib'),
    loadlibrary('invlib','invlib.h','alias','invlib'); % load the library
end

nullPtr = libpointer('voidPtr'); % empty pointer to void
fluxPtr = libpointer('doublePtr',flux); % pointer to one double
dlogfluxPtr = libpointer('doublePtr',dlogflux); % pointer to one double

result = calllib('invlib','ana_spec_inv',c,dc,Egrid,H,b,int_params,real_params,nullPtr,Eout,fluxPtr,dlogfluxPtr);
disp(sprintf('%s: ana_spec_inv result= %i',mfilename,result));

flux = fluxPtr.value;
dlogflux = dlogfluxPtr.value;

figure;
loglog(Eout,flux,'ro-',Eout,flux.*exp(-dlogflux*2),'ro--',Eout,flux.*exp(+dlogflux*2),'ro--');
set(gca,'xlim',[0.9,8.5]);
xlabel('Energy, MeV');
ylabel('Flux, #/cm^2/s/sr/MeV');

