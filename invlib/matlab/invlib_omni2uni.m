% demo of omni2uni, same as omni2uni_test.c

omniflux = 1.0E4; % typical value 1E+4 */
dlogomniflux = log(2)/2;
E = 300; % keV
BB0 = 4000.0/100.0; % B/B0
Lm = 6.6;

int_params = nan(5,1);
real_params = nan(3,1);

int_params(1+0) = 50; % NA - number of angular gridpoints */
int_params(1+1) = -1; % TEM-1 method */
int_params(1+2) = 0; % verbose = 0, not text output */
int_params(1+3) = 0; % root finder, default */
int_params(1+4) = 10e3; % maximumn # of iterations */

real_params(1+0) = E; % 300 keV */
real_params(1+1) = BB0; % B/B0 */
real_params(1+2) = Lm; % Lm */


if ~libisloaded('invlib'),
    loadlibrary('invlib','invlib.h','alias','invlib'); % load the library
end

nullPtr = libpointer('cstring'); % empty pointer to char *
unifluxPtr = libpointer('doublePtr',0); % pointer to one double
dlogunifluxPtr = libpointer('doublePtr',0); % pointer to one double

result = calllib('invlib','omni2uni',omniflux,dlogomniflux,int_params,real_params,nullPtr,unifluxPtr,dlogunifluxPtr);

uniflux = unifluxPtr.value;
dloguniflux = dlogunifluxPtr.value;

fprintf('---------------Manual lib call----------------\n');
fprintf('%s:inputs = omni=%.5e (dlog=%.5e) @ [%.1f keV, B/B0=%.2f, Lm=%.2f]\n',...
    mfilename,omniflux,dlogomniflux,E,BB0,Lm);
fprintf('%s:result = %i, uni=%.5e (dlog=%.5e)\n',mfilename,result,uniflux,dloguniflux);

fprintf('---------------Wrapper call----------------\n');

[uniflux,dloguniflux,result] = invlib('omni2uni',omniflux,dlogomniflux,'TEM1',...
    'keV',E,'B/B0',BB0,'Lm',Lm);
fprintf('%s:inputs = omni=%.5e (dlog=%.5e) @ [%.1f keV, B/B0=%.2f, Lm=%.2f]\n',...
    mfilename,omniflux,dlogomniflux,E,BB0,Lm);
fprintf('%s:result = %i, uni=%.5e (dlog=%.5e)\n',mfilename,result,uniflux,dloguniflux);

fprintf('---------------Multi call----------------\n');
s = (90:110)'/100;
[uniflux,dloguniflux,result] = invlib('omni2uni',omniflux*s,dlogomniflux,'TEM1',...
    'keV',E*s,'B/B0',BB0*s,'Lm',Lm*s);
for i = 1:length(s),
    fprintf('s=%g:result = %i, uni=%.5e (dlog=%.5e)\n',s(i),result(i),uniflux(i),dloguniflux(i));    
end
