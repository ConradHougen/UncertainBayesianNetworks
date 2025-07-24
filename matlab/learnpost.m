function aout = learnpost(nin)

Na = nin.Na;
Tx = nin.Tx;
Ny = nin.Ny;
Tyx = nin.Tyx;
Tynx = nin.Tynx;
nx = nin.nx;
nyx = nin.nyx;
nynx = nin.nynx;
ny = nin.ny;

ny1 = ny;
ny0 = Ny-ny1;

nx1 = nx;
nx0 = Tx-nx;

nyx1 = nyx;
nyx0 = Tyx-nyx1;
nynx1 = nynx;
nynx0 = Tynx-nynx1;

pex = nx1/(Tx+2);
peyx = nyx1/(Tyx+2);
peynx = nynx1/(Tynx+2);R = diag([pex*(1-pex)/Tx peyx*(1-peyx)/Na/pex peynx*(1-peynx)/Na/(1-pex)]);
    
nx1 = nx1+1;
nx0 = nx0+1;
nyx0 = nyx0+1;
nyx1 = nyx1+1;
nynx0 = nynx0+1;
nynx1 = nynx1+1;
    

e = [peyx-peynx pex 1-pex]';
pey = peyx*pex+peynx*(1-pex);
aout.ay_full = [pey 1-pey]*(pey*(1-pey)/(e'*R*e)-1);
    
v1 = (0:ny1);
c1 = ones(ny0+1,1)*(1./beta(v1+1,ny1-v1+1)/(ny1+1));
v1 = ones(ny0+1,1)*v1;
v1r = v1(:,end:-1:1);
    
v0 = (0:ny0)';
c0 = (1./beta(v0+1,ny0-v0+1)/(ny0+1))*ones(1,ny1+1);
v0 = v0*ones(1,ny1+1);
v0r = v0(end:-1:1,:);
    
denom = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)));
mx = sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
myx = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
mynx = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r+1,nynx0+v0r)))/denom;
m = [mx; myx; mynx];
r(1,1) = sum(sum(c1.*c0.*beta(nx1+v1+v0+2,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
r(2,2) = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1+2,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
r(3,3) = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r+2,nynx0+v0r)))/denom;
r(1,2) = sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
r(2,1) = r(1,2);
r(1,3) = sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r+1,nynx0+v0r)))/denom;
r(3,1) = r(1,3);
r(2,3) = sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r+1,nynx0+v0r)))/denom;
r(3,2) = r(2,3);
r = r-m*m';
    
my = (sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))+sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r+1).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r+1,nynx0+v0r))))/denom;
vy = (sum(sum(c1.*c0.*beta(nx1+v1+v0+2,nx0+v1r+v0r).*beta(nyx1+v1+2,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))+sum(sum(c1.*c0.*beta(nx1+v1+v0,nx0+v1r+v0r+2).*beta(nyx1+v1,nyx0+v0).*beta(nynx1+v1r+2,nynx0+v0r)))+2*sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r+1).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r+1,nynx0+v0r))))/denom;
vy = vy-my^2;
    
aout.ay2 = [my 1-my]*max(my*(1-my)/vy-1,2);
  
e = [myx-mynx mx 1-mx]';

my = myx*mx+mynx*(1-mx);

   
aout.ax = [mx 1-mx]*max(mx*(1-mx)/r(1,1)-1,2);
aout.ayx = [myx 1-myx]*max(myx*(1-myx)/r(2,2)-1,2);
aout.aynx = [mynx 1-mynx]*max(mynx*(1-mynx)/r(3,3)-1,2);
aout.ay = [my 1-my]*(my*(1-my)/(e'*r*e)-1);
    

a = myx*mx;
ea = [myx mx 0]';
va = ea'*r*ea;
b = mynx*(1-mx);
eb = [-mynx 0 1-mx]';
vb = eb'*r*eb;
vab = eb'*r*ea;
mxy = a/(a+b);
vxy = (b^2*va+a^2*vb-2*a*b*vab)/(a+b)^4;
aout.axy = [mxy 1-mxy]*max(mxy*(1-mxy)/vxy-1,2);
va =  ea'*diag(diag(R))*ea;
vb = eb'*diag(diag(R))*eb;
vab = eb'*diag(diag(R))*ea;
vxy = (b^2*va+a^2*vb-2*a*b*vab)/(a+b)^4;
aout.axy_ind = [mxy 1-mxy]*max(mxy*(1-mxy)/vxy-1,2);

   

a = (1-myx)*mx;
ea = [1-myx -mx 0]';
va = ea'*r*ea;
b = (1-mynx)*(1-mx);
eb = [-(1-mynx) 0 -(1-mx)]';
vb = eb'*r*eb;
vab = eb'*r*ea;
mxny = a/(a+b);
vxny = (b^2*va+a^2*vb-2*a*b*vab)/(a+b)^4;
aout.axny = [mxny 1-mxny]*max(mxny*(1-mxny)/vxny-1,2);
va =  ea'*diag(diag(R))*ea;
vb = eb'*diag(diag(R))*eb;
vab = eb'*diag(diag(R))*ea;
vxny = (b^2*va+a^2*vb-2*a*b*vab)/(a+b)^4;
aout.axny_ind = [mxny 1-mxny]*max(mxny*(1-mxny)/vxny-1,2);

v1_old = v1;
c1_old = c1;
v1r_old = v1r;
if (ny1>1),
       
    v1 = (0:ny1-1);
    c1 = ones(ny0+1,1)*(1./beta(v1+1,ny1-v1)/(ny1));
    v1 = ones(ny0+1,1)*v1;
    v1r = v1(:,end:-1:1);
    v0 = v0(:,1:end-1);
    c0 = c0(:,1:end-1);
    v0r = v0r(:,1:end-1);
    mxy = sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1+1,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom; 
    v1 = (0:ny1-2);
    c1 = ones(ny0+1,1)*(1./beta(v1+1,ny1-v1-1)/(ny1-1));
    v1 = ones(ny0+1,1)*v1;
    v1r = v1(:,end:-1:1);
    v0 = v0(:,1:end-1);
    c0 = c0(:,1:end-1);
    v0r = v0r(:,1:end-1);
    vxy = sum(sum(c1.*c0.*beta(nx1+v1+v0+2,nx0+v1r+v0r).*beta(nyx1+v1+2,nyx0+v0).*beta(nynx1+v1r,nynx0+v0r)))/denom;
    vxy = vxy-mxy^2;
    aout.axy2 = [mxy 1-mxy]*max(mxy*(1-mxy)/vxy-1,2);
else
    aout.axy2 = [NaN NaN];
end

if ny0>1,
    v1 = v1_old(1:end-1,:);
    c1 = c1_old(1:end-1,:);
    v1r = v1r_old(1:end-1,:);
    v0 = (0:ny0-1)';
    c0 = (1./beta(v0+1,ny0-v0)/(ny0))*ones(1,ny1+1);
    v0 = v0*ones(1,ny1+1);
    v0r = v0(end:-1:1,:);
    mxny = sum(sum(c1.*c0.*beta(nx1+v1+v0+1,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0+1).*beta(nynx1+v1r,nynx0+v0r)))/denom;
    v0 = (0:ny0-2)';
    c0 = (1./beta(v0+1,ny0-v0-1)/(ny0-1))*ones(1,ny1+1);
    v0 = v0*ones(1,ny1+1);
    v0r = v0(end:-1:1,:);
    v1 = v1_old(1:end-2,:);
    c1 = c1_old(1:end-2,:);
    v1r = v1r_old(1:end-2,:);
    vxny = sum(sum(c1.*c0.*beta(nx1+v1+v0+2,nx0+v1r+v0r).*beta(nyx1+v1,nyx0+v0+2).*beta(nynx1+v1r,nynx0+v0r)))/denom;
    vxny = vxny-mxny^2;
    aout.axny2 = [mxny 1-mxny]*max(mxny*(1-mxny)/vxny-1,2);
else
    aout.axny2 = [NaN NaN];
end
    








