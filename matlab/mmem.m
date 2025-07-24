function [p,cov,alpha] = mmem(circ,meta,val)

nparam = size(meta.varmap,2);
nnodes = size(meta.valuemap,2);

alpha = ones(nparam,1);

N = size(meta.valuemap,1);

M = ones(nparam,1)*meta.varmap;
maskfull = M==M';

%card = cell2mat({node(:).cardinality});
card = max(meta.valuemap,[],2);

cnt = 0;
mask = zeros(nparam,nparam);
for i=1:N
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        cnt = cnt+card(i);
    end
end

M1 = eye(nparam);
M2 = mask-M1;
M3 = maskfull-mask;

Ntrain = size(val,1);

S = mask*alpha;
p = alpha./S;

for i=1:Ntrain

    meta.p = p;
    obs = [(1:N)' val(i,:)'];
    obs = obs(val(i,:)>0,:);
    [~,m] = inferProb(circ,meta,obs);
    m2 = m(1:nparam).*m((1:nparam)+nnodes);

    p1 = p*ones(1,nparam);
    p2 = ((alpha+1)./(S+1))*ones(1,nparam);
    p22 = (alpha./(S+1))*ones(1,nparam);
    p3 = ((alpha+2)./(S+2))*ones(1,nparam);
    p32 = ((alpha+1)./(S+2))*ones(1,nparam);

    mn = ((p2.*M1+p22.*M2+p1.*M3)*m2)/m(nnodes);
    vv = ((p3.*p2.*M1+p32.*p22.*M2+p2.*p1.*M3)*m2)/m(nnodes);

    S = (mask*(mn.*(1-mn).*(mn-vv)))./(mask*(mn.*(1-mn).*(vv-mn.^2)));
    alpha = mn.*S;
    p = mn;

end


% EM PART

Nparam = length(meta.varmap);

map = meta.valuemap(:,1:Nparam);
card = max(map,[],2);

Nnodes = size(meta.valuemap,2);

nsee = sum(val>0,2);
valc = val(nsee==N,:);
valinc = val(nsee<N,:);

cs = [1; cumprod(card(1:end-1))];

theta0 = ones(Nparam,1);
for i=1:Nparam
    idx = max(meta.valuemap(:,i)-1,0)'*cs;
    vc = valc;
    vc(:,meta.valuemap(:,i)==0) = 0;
    theta0(i) = theta0(i)+sum(max(vc-1,0)*cs==idx);
end

cnt = 0;
cnt2 = 0;
mask = zeros(Nparam,Nparam);
for i=1:N
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        mask2(cnt2+1:cnt2+card(i)-1,cnt+1:cnt+card(i)) = [eye(card(i)-1) -ones(card(i)-1,1)];
        cnt = cnt+card(i);
        cnt2 = cnt2+card(i)-1;
    end
end

Nfree = size(mask2,1);

meta.p = p;
theta = theta0;
Nvalinc = size(valinc,1);

H = zeros(Nfree,Nfree);   
for i=1:Nvalinc
    obs = [(1:N)' valinc(i,:)'];
    obs = obs(obs(:,2)>0,:);
    [~,m] = inferProb(circ,meta,obs);
    theta = theta+m(1:Nparam).*m(Nnodes+(1:Nparam))/m(Nnodes);
    v = mask2*((m(1:Nparam)>0).*m(Nnodes+(1:Nparam)));
    H = H+v*v'/m(Nnodes)^2;
end

Z = mask*theta;
meta.p = theta./Z;

for i=1:size(valc,1)
    v = mask2*((sum((valc(i,:)'*ones(1,Nparam)).*(meta.valuemap(:,1:Nparam)>0)==meta.valuemap(:,1:Nparam),1)==N)'./meta.p);
    H = H + v*v';
end

p = meta.p;
cov = inv(eye(Nfree)+H);
cov = mask2'*cov*mask2;

