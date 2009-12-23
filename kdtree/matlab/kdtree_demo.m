% demonstrate kdtree nearest neighbors

%% small 2-d test
disp('starting 2-d test');
X = rand(1e4,2);
X0 = rand(1,size(X,2));
tic;
tree = kdtree('build',X);
toc;
tic;
[index,R2,XNN] = kdtree('kNN',X,tree,X0,30);
toc;
r = sqrt(R2(end));
t = (0:100)*2*pi/100;
figure;
plot(X(:,1),X(:,2),'k.',...
    X0(1),X0(2),'bs',...
    X0(1)+cos(t)*r,X0(2)+sin(t)*r,'b-',...
    X(index,1),X(index,2),'ro',...
    XNN(:,1),XNN(:,2),'rs');
axis([0 1 0 1]);
title('2-D (kNN)');

% small 2-d test of lib
disp('calling lib (2d)');
X0 = rand(5,size(X,2));
tic;
[index,R2,XNN] = kdtree('multikNN',X,tree,X0,20);
toc;
t = (0:100)*2*pi/100;
figure;
plot(X(:,1),X(:,2),'k.');
hold on;
r2 = nan(size(R2));
for iX0 = 1:size(X0,1),
    r = sqrt(R2(end,iX0));
    h = plot(X0(iX0,1),X0(iX0,2),'*',...
        X0(iX0,1)+cos(t)*r,X0(iX0,2)+sin(t)*r,'-',...
        X(index(:,iX0),1),X(index(:,iX0),2),'o',...
        XNN(:,iX0,1),XNN(:,iX0,2),'s');
    set(h,'color',getcolor(iX0,size(X0,1)));
    d = sum((repmat(X0(iX0,:),size(X,1),1)-X).^2,2);
    [dsort,isort] = sort(d);
    isort = isort(1:size(index,1));
    plot(X(isort,1),X(isort,2),'h','color',getcolor(iX0,size(X0,1)),'markersize',10);
end
axis([0 1 0 1]);
title('2-D (multikNN)');

disp('checking for errors');

figure;
sumerrors = 0;
for iX0 = 1:size(X0,1),
    d = sum((repmat(X0(iX0,:),size(X,1),1)-X).^2,2);
    [dsort,isort] = sort(d);
    plot(index(:,iX0),isort(1:size(index,1)),'.','color',getcolor(iX0,size(X0,1)));
    hold on;
    errors = sum(abs(double(index(:,iX0))-isort(1:length(index))));
    sumerrors = sumerrors + errors;
    if errors,
        disp(iX0);
        testinfo.iX0 = iX0;
        testinfo.X = X;
        testinfo.X0 = X0;
        testinfo.tree = tree;
        testinfo.k = size(index,1);
        testinfo.dsort = dsort;
        testinfo.isort = isort;
        testinfo.d = d;
        testinfo.index = index;
        testinfo.R2 = R2;
        testinfo.XNN = XNN;
        [test2.index,test2.R2,test2.XNN] = kdtree('multikNN',testinfo.X,testinfo.tree,testinfo.X0,testinfo.k);
        [testinfo.one.index,testinfo.one.R2,testinfo.one.XNN] = kdtree('kNN',testinfo.X,testinfo.tree,testinfo.X0(testinfo.iX0,:),testinfo.k);
        plot(testinfo.index(:,testinfo.iX0),testinfo.isort(1:testinfo.k),'o','color',getcolor(iX0,size(X0,1)))
        plot(test2.index(:,testinfo.iX0),testinfo.isort(1:testinfo.k),'s','color',getcolor(iX0,size(X0,1)))
        plot(testinfo.one.index,testinfo.isort(1:testinfo.k),'+','color',getcolor(iX0,size(X0,1)))
        testinfo.old.tree = kdtree_old('build',testinfo.X);
        [testinfo.old.index,testinfo.old.R2,testinfo.old.XNN] = kdtree_old('kNN',testinfo.old.tree,testinfo.X0(testinfo.iX0,:),testinfo.k);
    end
end
title(sprintf('%d-D, index errors: %d',size(X0,2),sumerrors));


%%
% large 3-d test
disp('large, 3-d test');
k = 30;
X = rand(1e6,3);
tic;
tree = kdtree('build',X);
toc;
disp(sprintf('Finding nearest neighbors for %d points',size(X0,1)));
X0 = rand(1,size(X,2));
tic;
[index,R2,XNN] = kdtree('kNN',X,tree,X0,k);
toc;
d = sum((repmat(X0,size(X,1),1)-X).^2,2);
[dsort,isort] = sort(d);
figure;
plot(index,isort(1:length(index)),'.');

title(sprintf('3-D, index errors: %d',sum(abs(double(index)-isort(1:length(index))))));

disp('large, 3-d test of lib');
% large 3-d test of lib
X0 = [X0;rand(10000-1,size(X,2))];
disp(sprintf('Finding nearest neighbors for %d points',size(X0,1)));
tic;
[index,R2,XNN] = kdtree('multikNN',X,tree,X0,k);
toc;
figure;
sumerrors = 0;
for iX0 = 1:min(size(X0,1),20),
    d = sum((repmat(X0(iX0,:),size(X,1),1)-X).^2,2);
    [dsort,isort] = sort(d);
    plot(index(:,iX0),isort(1:size(index,1)),'.','color',getcolor(iX0,size(X0,1)));
    hold on;
    sumerrors = sumerrors+sum(abs(double(index(:,iX0))-isort(1:k)));
end
title(sprintf('3-D, index errors: %d',sumerrors));


disp('done');
