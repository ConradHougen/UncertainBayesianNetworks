function [p,cov,alpha] = mm_lbp(meta,node,val)

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
for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
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

for i=1:Ntrain,

    meta.p = p;
    %%%
    %[~,m] = inferProbfast(meta,val(i,:));
    node = prob_infer_loopy(node,val(i,:));
    % Concatenate all the nparam means
    m = zeros(np
    for j = 1:length(node)
        
    end
    
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

cov = (diag(p./(S+1))-(p./(S+1))*p').*mask;