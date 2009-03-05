% demo of wide2uni, same as omni2uni_test.c, but computes PAgrid, H, and
% sets ialpha0

int_params = nan(5,1);
real_params = nan(3,1);

NA = 50; % NA - number of angular gridpoints */

int_params(1+0) = NA; % NA - number of angular gridpoints */
int_params(1+1) = -1; % TEM-1 method */
int_params(1+2) = 0; % verbose = 0, not text output */
int_params(1+3) = 3; % minimizer, 0=BFGS, 3=NM */
int_params(1+4) = 1000; % maximumn # of iterations */

real_params(1+0) = 300.0; % 300 keV */
real_params(1+1) = 40000.0/100.0; % B/B0 */
real_params(1+2) = 6.6; % Lm */

wideflux = 1.0E4; % typical value 1E+4 */
dlogwideflux = log(2)/2;

% generate PAgrid, H, and ialpha0, just as in C omni2uni
da = 90/NA;
PAgrid = (1:NA)*da;
H = sin(PAgrid*pi/180)*da*pi/180;
H(end) = H(end)/2;
ialpha0 = length(PAgrid)-1; % locally mirroring

if ~libisloaded('invlib'),
    loadlibrary('invlib','invlib.h','alias','invlib'); % load the library
end

nullPtr = libpointer('cstring'); % empty pointer to char *
unifluxPtr = libpointer('doublePtr',0); % pointer to one double
dlogunifluxPtr = libpointer('doublePtr',0); % pointer to one double

result = calllib('invlib','wide2uni',wideflux,dlogwideflux,PAgrid,H,int_params,real_params,nullPtr,ialpha0,unifluxPtr,dlogunifluxPtr);

uniflux = unifluxPtr.value;
dloguniflux = dlogunifluxPtr.value;

disp(sprintf('%s:inputs = wide=%.5e (dlog=%.5e) @ [%.1f keV, B/B0=%.2f, Lm=%.2f]',...
    mfilename,wideflux,dlogwideflux, real_params(1+0),real_params(1+1),real_params(1+2)));
disp(sprintf('%s:result = %i, uni=%.5e (dlog=%.5e)',mfilename,result,uniflux,dloguniflux));
