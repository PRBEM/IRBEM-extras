% demo ana_spec_inv - just like specinv_test.c

infile = 'specinv_test.in1';

if ~exist(infile,'file'),
    error('Unable to locate "%s". Please copy it into matlab path or edit "%s.m"',infile,mfilename);
end

outfile = [which(infile),'_matlab_log'];

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

H = reshape(H,NE,NC)';

% don't need Eout

% load cress sample mean/cov data from CRRES MEA+HEEF combined
% for L~4-5

% Nq_max is the maximum number of bases you expect to use,
% Nq_max <= length(crres.MeV)
crres = load('crres_demo_data.mat');
Nq_max = 10;
[V,D] = eigs(crres.covX,Nq_max,'la',struct('issym',1,'disp',0));

% assume flux is zero outside crres energy range
iEgrid = find((Egrid>=crres.MeV(1)) & (Egrid<=crres.MeV(end)));
H = H(:,iEgrid);
Egrid = Egrid(iEgrid);
NE = length(Egrid);
% in practice, we might use AE8 to extend the energy range above/below the
% energy range of CRRES MEA+HEEF

% now, interpolate the mean_log_flux and the V's (defined on crres.MeV) onto the Egrid
mean_log_flux = interp1(log(crres.MeV),crres.meanX,log(Egrid),'linear');
basis_vectors = interp1(log(crres.MeV),V,log(Egrid),'linear');
basis_variance = diag(D);

figure;
subplot(2,1,1);
loglog(Egrid,exp(mean_log_flux),'k-','linew',2);
xlabel('MeV');
ylabel('#/cm^2/s/sr/MeV');
legend('<log flux>');
title('CRRES MEA/HEEF L~4-5 PC Model');
grid on;
subplot(2,2,3);
semilogx(Egrid,basis_vectors(:,1:3),'linew',2);
xlabel('MeV');
ylabel('Basis Vectors 1-3');
grid on;
subplot(2,2,4);
semilogy(basis_variance,'o-','linew',2);
xlabel('PC #');
ylabel('Variance in log flux');
grid on;

% now do pc fit
NQs = [3,5,Nq_max];
fits = cell(size(NQs));
for iNq = 1:length(NQs),
    Nq = NQs(iNq);
    % a reduced number of bases to use in the inversion
    fits{iNq} = invlib('pc_spec_inv',c(:)',dc(:)',Egrid,H,1,b(:)',mean_log_flux,basis_vectors,basis_variance,...
        'dE_mode','G=GdE','outfile',[outfile,'4'],'num_bases',Nq);
end

% for reference, the analytical fit
asi_fit = invlib('ana_spec_inv',c(:)',dc(:)',Egrid,H,1,b(:)',Egrid,'pl','exp','dE_mode','G=GdE','outfile',[outfile,'3']);

figure;
h1 = loglog(Egrid,asi_fit.flux,'ko-',Egrid,asi_fit.flux.*exp(-asi_fit.dlogflux*2),'ko--',Egrid,asi_fit.flux.*exp(+asi_fit.dlogflux*2),'ko--');
hold on;
h2 = loglog(Egrid,exp(mean_log_flux),'k-','linew',2);
leg = {'Analytical Fit','<log flux> from PC model'};
h = [h1(1),h2(1)];
for iNq = 1:length(NQs),
    fit = fits{iNq};
    h3 = loglog(Egrid,fit.flux,'.-',Egrid,fit.flux.*exp(-fit.dlogflux*2),'.--',Egrid,fit.flux.*exp(+fit.dlogflux*2),'.--');
    set(h3,'color',getcolor(iNq,length(NQs)));
    h(end+1) = h3(1);
    leg{end+1} = sprintf('%d-PC fit',length(fit.q));
end

%set(gca,'xlim',[0.9,8.5]);
xlabel('Energy, MeV');
ylabel('Flux, #/cm^2/s/sr/MeV');
title('ICO fit example');
legend(h,leg{:});

% % this code chunk tests calculation of dlogflux which, tragically
% % is NaN sometimes when the nonlinear optimization fails (the hessian is
% % not positive definite at the stopping point).
%
% flux = fit.flux';
% lambda = H*flux;
% logy = log(c(:));
% loglam = log(lambda);
% sigma = dc(:);
% z = (logy-loglam)./sigma;
% pen_ell = sum(0.5*z.^2);
% pen_grad = -z./sigma./lambda;
% pen_hess = (1+logy-loglam)./(sigma.*lambda).^2;
% dfdq = diag(sparse(fit.flux))*basis_vectors(:,1:Nq);
%
% Nq = length(fit.q);
% d2elldq2 = zeros(Nq);
% for m = 1:Nq,
%     for mp = 1:Nq,
%         for i = 1:NC,
%             for j = 1:NE,
%                 for jp = 1:NE,
%                     d2elldq2(m,mp) = d2elldq2(m,mp)+pen_hess(i)*H(i,j)*dfdq(j,m)*H(i,jp)*dfdq(jp,mp);
%                 end
%             end
%         end
%         for i = 1:NC,
%             for j = 1:NE,
%                 d2elldq2(m,mp) = d2elldq2(m,mp)+pen_grad(i)*H(i,j)*(dfdq(j,m)*basis_vectors(j,mp));
%             end
%         end
%     end
% end
% d2elldq2 = d2elldq2+diag(1./basis_variance(1:Nq));
%
% g = basis_vectors(:,1:Nq);
% varlogflux = g*inv(d2elldq2)*g';
% dlogflux = sqrt(diag(varlogflux));
% max(abs(dlogflux-fit.dlogflux'))
%
% disp([any(~isfinite(fit.dlogflux)) any(~isfinite(dlogflux))]);
%
% loglog(Egrid,fit.flux'.*exp(-dlogflux*2),'gs--',Egrid,fit.flux'.*exp(+dlogflux*2),'gs--');
%
% figure;
% plot(dlogflux,fit.dlogflux,'.');
