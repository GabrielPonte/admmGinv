m = 2000;
n = 1000;
r = 250;
dens = 0.25;
M = 2;

rho=(1/M)^(2/(r+1));

rc=M*rho.^(1:r);

A= full(sprand(m,n,dens,rc));

%[C,time] = greedy_light(A,m,n,r,[],1);
%[C2,time_new] = qrGL(A,m,n,r,1);
tic
[~,~,C3] = qr(A,"econ","vector");
C3 = C3(1:r);
time_qr = toc;
[time_new,time_qr]
[min(svd(A(:,C2))),min(svd(A(:,C3)))]