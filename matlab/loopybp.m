function [p,cov] = loopybp(meta,val)

nparam = length(meta.varmap);
[N,nnodes] = size(meta.valuemap);
Ntrain = size(val,1);

alpha = ones(nparam,1);

card = max(meta.valuemap,[],2);

M = ones(nparam,1)*meta.varmap;
maskfull = M==M';

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

S = mask*alpha;
p = alpha./S;

for i=1:Ntrain
    meta.p = p;
    obs = [(1:N)' val(i,:)'];
    obs = obs(val(i,:)>0,:);
    
    
    
end

end