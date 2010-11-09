% debug invlib errors reported by Steve Morley

% Hi Paul,
% The response functions (geometric factors as a function of energy) are 
% given in the file bdd2r_ns41_electrons_current.txt. The files should be 
% intelligible, but just in case the first column is the incident energy 
% [keV] and each column represents the calculated response [cm^2 sr] to 
% that energy for each channel. 
respfile = [finddrive('irbem'),filesep,'extras',filesep,'invlib',filesep,'debug',filesep,'bdd2r_ns41_electrons_current.txt'];
outfile= [finddrive('irbem'),filesep,'extras',filesep,'invlib',filesep,'debug',filesep,'test.out'];

% 9 colums: keV then 8 channels of cm^s sr
m0c2 = 511; % electron rest mass, keV

fid = fopen(respfile,'r');
C = textscan(fid,repmat('%f',1,9),'headerlines',13);
fclose(fid);
keV = C{1};
epsG = cat(2,C{2:end});
%The integration interval is 240s and this 
% has not been taken into account in this file - I do that in my interface.
dt = 240;
% Also, the output energy grid is defined (when I'm running the code, for 
% now at least) to just be the set of energies at which the response 
% functions are defined.
% 
% I'll give you two sets of counts and their estimated backgrounds. I've 
% already multiplied through by the 240 second accumulation interval so 
% they're counts/accumulation rather than a rate. In each case I've set dy 
% (relative error) as 0.3466 for each channel. I.e. dy = [0.3466, 0.3466, 
% 0.3466,
% 0.3466, 0.3466, 0.3466, 0.3466, 0.3466]
dy = repmat(0.3466,8,1);
% 
% If you need more info to attempt to reproduce what I'm doing, let me know.
% 
Case = 5;
switch(Case),
    case 1,
        % Case 1 (ns41 at L~9):
        Counts = [121659,  11173, 450317, 868, 765, 678, 769, 734];
        Background = [384.072, 408.192, 379.128, 390.624, 406.872, 413.424, 420.144, 391.248];
    case 2,
        % Case 2 (ns41 at L~4.3):
        Counts =[ 999283, 204579, 34407, 28297, 14745, 12804, 11814, 11020];
        Background=[ 384.192,  408.36 ,  378.984,  390.768,  406.896,  413.568, 420.336,  391.416];
    case 3, % fake power law
        Background=[ 384.192,  408.36 ,  378.984,  390.768,  406.896,  413.568, 420.336,  391.416];
        %Background= zeros(1,8)+rand(1,8);
        Counts = Background + dt*trapz(keV,epsG.*repmat(1e2*(keV/100).^-2.5,1,8));
        fprintf('%10f, %10f, %10f\n',[Background', Counts', Counts'-Background']')
    case 4, % fake channels and fake power law
        Background=[ 384.192,  408.36 ,  378.984,  390.768,  406.896,  413.568, 420.336,  391.416];
        epsG = double([ keV > 100 , keV > 200 , keV > 400 , keV > 800 , keV > 1600 , keV > 3200, keV > 6400, keV > 12800]);
        % Counts = Background + dt*trapz(keV,epsG.*repmat(1e4*(keV/100).^-2,1,8));
        Counts = Background + dt*trapz(keV,epsG.*repmat(exp(10-3*log(keV)),1,8));
        fprintf('%10f, %10f, %10f\n',[Background', Counts', Counts'-Background']')
    case 5, % fake rm2
        Background=[ 384.192,  408.36 ,  378.984,  390.768,  406.896,  413.568, 420.336,  391.416];
        q = [2 -1/100 2 -1/300];
        rm2 = @(q,E)E.*(1+E/2/m0c2).*exp(q(1)+E*q(2))+E.*(1+E/2/m0c2).*exp(q(3)+E*q(4));
        Counts = Background + dt*trapz(keV,epsG.*repmat(rm2(q,keV),1,8));
        fprintf('%10e, %10e, %10e\n',[Background', Counts', Counts'-Background']')
end
% 
% These break at my end unless I use the PLE function and a minimizer 
% other than the BFGS.
% 
% A commented Python script is also provided that shouldn't need anything 
% beyond Python and Numpy installed. If you have Matplotlib then that'll 
% let you use the plot() method to quickly visualize the result. Oh, 
% you'll need to install pyinvlib from the IRBEM svn too.
% That's as simple as going to the invlib directory within your local copy 
% of the repo and typing:
% sudo python setup.py install
% 
% Then from the directory that the Python script is in, type:
% python ns41_test.py
% and everything should just work (assuming you have numpy (and 
% matplotlib) and pyinvlib installed -- I've not tested pyinvlib under 
% cygwin, so this may require a linux box (or VM) to get going. Let me 
% know how it works under cygwin (if you try that).
% If you have ipython, use that and then type:
% run -i ns41_test.py
% so you can play with things after running the script.
% 
% 
% Cheers,
% Steve
% 
% PS. Channel 3 is known to be intermittently dodgy at this time, so I've 
% tried excluding that, but that still hasn't made things work for me.


ii = [1:8]; % leave out first channel
fit = invlib('ana_spec_inv',Counts(ii),dy(ii),keV,epsG(:,ii)',dt,Background(ii),keV,'trapz','RM2','outfile',outfile,'BFGS','rest_energy',m0c2);
subplot(2,1,1);
loglog(keV,fit.flux);
subplot(2,1,2);
ax = [min(Counts(ii)) max(Counts(ii))];
loglog(ax,ax,'k-',Counts(ii),fit.lambda,'bo',Background(ii),fit.lambda,'rs');

