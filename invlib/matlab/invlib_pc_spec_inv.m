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

% don't need Eout
asi_fit = invlib('ana_spec_inv',c(:)',dc(:)',Egrid,reshape(H,NE,NC)',1,b(:)',Egrid,'pl','exp','dE_mode','G=GdE','outfile',[outfile,'3']);
V = [ones(NE,1) (1:NE)',((0.5:NE)'-(NE/2)).^2]; % flat, linear, and quadratic PCs
V(:,1) = V(:,1)/norm(V(:,1));
V(:,2) = V(:,2)-mean(V(:,2));
V(:,2) = V(:,2)/norm(V(:,2));
V(:,3) = V(:,3)-mean(V(:,3));
V(:,3) = V(:,3)/norm(V(:,3));
D = diag([3;2;1]); % dummy values for the 

basis_vectors = V;
basis_variance = diag(D);
% fit = invlib('pc_spec_inv',y,dy,Egrid,G,dt,b,mean_log_flux,basis_vectors,basis_variance...)

Nq = 3;

% force the mean log flux to a random deviation from the analytical
% solution, plus some noise
mean_log_flux = log(asi_fit.flux)'+V*(randn(size(V,2),1).*sqrt(basis_variance))+randn(size(asi_fit.flux'))/10; 

fit = invlib('pc_spec_inv',c(:)',dc(:)',Egrid,reshape(H,NE,NC)',1,b(:)',mean_log_flux,basis_vectors,basis_variance,...
    'dE_mode','G=GdE','outfile',[outfile,'4'],'num_bases',Nq);

figure;
h1 = loglog(Egrid,asi_fit.flux,'b.-',Egrid,asi_fit.flux.*exp(-asi_fit.dlogflux*2),'b.--',Egrid,asi_fit.flux.*exp(+asi_fit.dlogflux*2),'b.-');
hold on;
h2 = loglog(Egrid,exp(mean_log_flux),'k-');
h3 = loglog(Egrid,fit.flux,'ro-',Egrid,fit.flux.*exp(-fit.dlogflux*2),'ro--',Egrid,fit.flux.*exp(+fit.dlogflux*2),'ro--');
set(gca,'xlim',[0.9,8.5]);
xlabel('Energy, MeV');
ylabel('Flux, #/cm^2/s/sr/MeV');
legend([h1(1),h2(1),h3(1)],'Analytical Fit','Assumed q=0 PC model','PC fit');

% stop
% this code chunk tests
% 
% H = reshape(H,length(Egrid),length(c))';
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
% Ny = length(c);
% d2elldq2 = zeros(Nq);
% for m = 1:Nq,
%     for mp = 1:Nq,
%         for i = 1:Ny,
%             for j = 1:NE,
%                 for jp = 1:NE,
%                     d2elldq2(m,mp) = d2elldq2(m,mp)+pen_hess(i)*H(i,j)*dfdq(j,m)*H(i,jp)*dfdq(jp,mp);
%                 end
%             end
%         end
%         for i = 1:Ny,
%             for j = 1:NE,
%                 d2elldq2(m,mp) = d2elldq2(m,mp)+pen_grad(i)*H(i,j)*(dfdq(j,m)*basis_vectors(j,mp));
%             end
%         end
%     end
% end
% 
% g = basis_vectors(:,1:Nq);
% varlogflux = g*(d2elldq2\g');
% varlogflux2 = g*inv(d2elldq2)*g';
% dlogflux = sqrt(diag(varlogflux));
% max(abs(dlogflux-fit.dlogflux'))
% 
% disp([any(~isfinite(fit.dlogflux)) any(~isfinite(dlogflux))]);
% 
% loglog(Egrid,fit.flux'.*exp(-dlogflux*2),'gs--',Egrid,fit.flux'.*exp(+dlogflux*2),'gs--');
% 
% figure;
% plot(dlogflux,fit.dlogflux,'.');
