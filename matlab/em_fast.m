function [pout,cov,mask] = em_fast(circ,meta,val,p0)

Nparam = length(meta.varmap);

map = meta.valuemap(:,1:Nparam);
card = max(map,[],2);

N = size(map,1); 
Nnodes = size(meta.valuemap,2);

nsee = sum(val>0,2);
valc = val(nsee==N,:);
valinc = val(nsee<N,:);

cs = [1; cumprod(card(1:end-1))];



theta0 = ones(Nparam,1);
for i=1:Nparam,
    idx = max(meta.valuemap(:,i)-1,0)'*cs;
    vc = valc;
    vc(:,meta.valuemap(:,i)==0) = 0;
    theta0(i) = theta0(i)+sum(max(vc-1,0)*cs==idx);
end



    

cnt = 0;
cnt2 = 0;
mask = zeros(Nparam,Nparam);
for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        mask2(cnt2+1:cnt2+card(i)-1,cnt+1:cnt+card(i)) = [eye(card(i)-1) -ones(card(i)-1,1)];
        cnt = cnt+card(i);
        cnt2 = cnt2+card(i)-1;
    end
end

Nfree = size(mask2,1);

Ntrain = size(valinc,1);




% p = ones(Nparam,1);
% k = mask*p;
% meta.p = p./k;

if nargin<4,
    meta.p = theta0./(mask*theta0);
else
    meta.p = p0;
end

go = 1;
cnt = 1;
err = 1;
while go,
    pold = meta.p;
    %hd = waitbar(0);
    theta = theta0;
    H = zeros(Nfree,Nfree);   
    for i=1:Ntrain,
        %waitbar(i/Ntrain,hd,sprintf('err = %f',err));
        obs = [(1:N)' valinc(i,:)'];
        obs = obs(obs(:,2)>0,:);
        [~,m] = inferProb(circ,meta,obs);
        theta = theta+m(1:Nparam).*m(Nnodes+(1:Nparam))/m(Nnodes);
        v = mask2*((m(1:Nparam)>0).*m(Nnodes+(1:Nparam)));
        H = H+v*v'/m(Nnodes)^2;
    end
    %delete(hd);
    Z = mask*theta;
    meta.p = theta./Z;
    cnt = cnt+1;
    err = sqrt(mean((meta.p-pold).^2));
    if (err<1e-3)||(cnt>=50),
        go = 0;
    end
end

for i=1:size(valc,1),
    v = mask2*((sum((valc(i,:)'*ones(1,Nparam)).*(meta.valuemap(:,1:Nparam)>0)==meta.valuemap(:,1:Nparam),1)==N)'./meta.p);
    H = H + v*v';
end

pout = meta.p;
cov = inv(eye(Nfree)+H);
cov = mask2'*cov*mask2;
    



        
    