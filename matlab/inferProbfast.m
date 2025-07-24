function [p,m] = inferProbfast(meta,val)

N = size(meta.valuemap,1);
Ng = length(meta.varmap);
card = max(meta.valuemap,[],2);

nnodes = size(meta.valuemap,2);

Inpt = zeros(Ng,size(val,1));
for i=1:N,
    Inpt = Inpt + meta.Li{i}(:,val(:,i)'+1);
end

m = zeros(2*nnodes,size(val,1));
m(1:Ng,:) = Inpt.*(meta.p*ones(1,size(val,1)));
cnt2 = 0;
for i=1:meta.Nlvl,
    mi = m(meta.liptf{i},:);
    if cnt2,
        mi(mi==0) = eps;  %avoid zero times inf
        mi = log(mi);
        m(meta.levelf{i},:) = exp(meta.Lf{i}*mi);
    else
        m(meta.levelf{i},:) = meta.Lf{i}*mi;
    end  
    cnt2 = xor(cnt2,1);
end

m(2*nnodes,:) = 1;
cnt2 = mod(meta.Nlvl+1,2);
for i=meta.Nlvl:-1:1,
    mi = m(meta.liptb{i},:);
    if cnt2,
        mi(mi==0) = eps;  %avoid zero times inf
        mi = log(mi);
        m(meta.levelb{i},:) = m(meta.levelb{i},:)+meta.Lb2{i}*exp(meta.Lb{i}*mi); 
    else
        m(meta.levelb{i},:) = m(meta.levelb{i},:)+meta.Lb{i}*mi;
    end
    cnt2 = xor(cnt2,1);
end

denom = m(nnodes,:);
p = NaN*zeros(N,max(card),size(val,1));
for i=1:N,
    for j=1:card(i),
        loc = find(meta.varmap==i);
        loc2 = loc(meta.valuemap(i,loc)==j);
        p(i,j,:) = sum(m(loc2,:).*m((loc2)+nnodes,:),1)./denom;
    end
end