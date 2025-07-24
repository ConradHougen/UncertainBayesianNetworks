function [p,v,mm] = inferMonte(net,obs,Nmonte)

if nargin<3,
    Nmonte = 100;
end

[circ,meta] = net2circ_monte(net);

alpha = meta.alpha;
alphaloc = meta.alphaloc;

n = size(meta.valuemap,1);
card = max(meta.valuemap,[],2);
mm = zeros(size(meta.valuemap,2)*2,Nmonte);

pm = zeros(n,max(card),Nmonte);
hd = waitbar(0);
for i=1:Nmonte,
    waitbar(i/Nmonte,hd);
    for l=alphaloc,
        pp = gamrnd(alpha{l},1);
        pp = pp/sum(pp);
        ln = length(alpha{l});
        meta.m(l:l+ln-1) = pp';
    end
    meta.p = meta.m;
    [pm(:,:,i), mm(:,i)] = inferProb(circ,meta,obs);
end
v = mean(pm.^2,3);
p = mean(pm,3);
v = v-p.^2;

% mv = mean(mm.^2,2);
% mm = mean(mm,2);
% mv = mv-mm.^2;

delete(hd);
    
    