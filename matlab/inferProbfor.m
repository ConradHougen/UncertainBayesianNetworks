function [p,m] = inferProbfor(meta,val)

N = size(meta.valuemap,1);
Ng = length(meta.varmap);
card = max(meta.valuemap,[],2);

nnodes = size(meta.valuemap,2);

Inpt = zeros(Ng,size(val,1));
for i=1:N,
    Inpt = Inpt + meta.Li{i}(:,val(:,i)'+1);
end

m = zeros(nnodes,size(val,1));
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

p = exp(sum(log(m(nnodes,:))));