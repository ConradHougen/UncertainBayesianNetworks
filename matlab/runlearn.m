%initally use the ground turth px instead of the estimated pex in EM, which
%resulted in a better estimated as predicted by the CRLB estimate of the
%variance.  The non-biased version of the EM means the variance is smaller
%than estimation variance.  Using the biased EM, seems to lower the
%variance of estimator to better match the CRLB based variance. 

tol = 1e-8;

nin.tol = tol;

Nruns = 1000;

Na = 10;
Nx = 0;
Ny = 100;

N = Na+Nx+Ny;

ax = zeros(Nruns,2);
ay = zeros(Nruns,2);
ay2 = zeros(Nruns,2);
ayx = zeros(Nruns,2);
aynx = zeros(Nruns,2);
axy = zeros(Nruns,2);
axny = zeros(Nruns,2);
axy2 = zeros(Nruns,2);
axny2 = zeros(Nruns,2);

ax_cr = zeros(Nruns,2);
ay_cr = zeros(Nruns,2);
ayx_cr = zeros(Nruns,2);
aynx_cr = zeros(Nruns,2);
axy_cr = zeros(Nruns,2);
axny_cr = zeros(Nruns,2);
axy_ind = zeros(Nruns,2);
axny_ind = zeros(Nruns,2);
axy_indcr = zeros(Nruns,2);
axny_indcr = zeros(Nruns,2);


pgtx = zeros(Nruns,1);
pgty = zeros(Nruns,1);
pgtyx = zeros(Nruns,1);
pgtynx = zeros(Nruns,1);
pgtxy = zeros(Nruns,1);
pgtxny = zeros(Nruns,1);

nsv = zeros(Nruns,4);

ax_full = zeros(Nruns,2);
ayx_full = zeros(Nruns,2);
aynx_full = zeros(Nruns,2);
ay_full = zeros(Nruns,2);

hd = waitbar(0);

% pyx = 1;
% pynx = 0;

for i=1:Nruns,
    
    waitbar(i/Nruns,hd);

    px = rand(1);
    pyx = rand(1);
    %pynx = pyx;
    pynx = rand(1);

    x = rand(N,1)<=px;
    y1 = rand(N,1)<=pyx;
    y2 = rand(N,2)<=pynx;
    y = y1;
    y(x==0) = y2(x==0);

    Tx = Na+Nx;
    nx = sum(x(1:Tx));
    nyx = sum(y(x(1:Na)==1));
    Tyx = sum(x(1:Na)==1);
    Tynx = Na-Tyx;
    nynx = sum(y(x(1:Na)==0));

    ny = sum(y(Tx+1:end));
    Ty = Ny;
    
    nin.Na = Na;
    nin.Tx = Tx;
    nin.Ty = Ty;
    nin.Ny = Ny;
    nin.Tyx = Tyx;
    nin.Tynx = Tynx;
    nin.nx = nx;
    nin.nyx = nyx;
    nin.nynx = nynx;
    nin.ny = ny;
    
    nsv(i,:) = [nx nyx nynx ny];
    
    aout = learnpost(nin);
    acr = learnb2(nin);

    
    ax_full(i,:) = [nx+1 Tx-nx+1];
    ayx_full(i,:) = [nyx+1 Tyx-nyx+1];
    aynx_full(i,:) = [nynx+1 Tynx-nynx+1];
    
    ax(i,:) = aout.ax;
    ayx(i,:) = aout.ayx;
    aynx(i,:) = aout.aynx;
    ay(i,:) = aout.ay;
    ay2(i,:) = aout.ay2;
    axy(i,:) = aout.axy;
    axny(i,:) = aout.axny;
    axy2(i,:) = aout.axy2;
    axny2(i,:) = aout.axny2;
    axy_ind(i,:) = aout.axy_ind;
    axny_ind(i,:) = aout.axny_ind;
    
    ax_cr(i,:) = acr.ax;
    ayx_cr(i,:) = acr.ayx;
    aynx_cr(i,:) = acr.aynx;
    ay_cr(i,:) = acr.ay;
    axy_cr(i,:) = acr.axy;
    axny_cr(i,:) = acr.axny;
    axy_indcr(i,:) = acr.axy_ind;
    axny_indcr(i,:) = acr.axny_ind;
    
    
    
    pgtx(i) = px;
    %pgtx(i) = pex;
    pgtyx(i) = pyx;
    pgtynx(i) = pynx;
    pgty(i) = pyx*px+pynx*(1-px);
    pgtxy(i) = 1/(1+pynx*(1-px)/pyx/px);
    pgtxny(i) = 1/(1+(1-pynx)*(1-px)/(1-pyx)/px);
    
end

delete(hd);
    
    








