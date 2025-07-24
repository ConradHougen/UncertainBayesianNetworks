function [p,v,m,cov] = inferSOProbfast(meta,val)

N = size(meta.valuemap,1);
Ng = length(meta.varmap);
card = max(meta.valuemap,[],2);

nnodes = size(meta.valuemap,2);

Inpt = zeros(Ng,1);
for i=1:N,
    Inpt = Inpt + meta.Li{i}(:,val(i)+1);
end

m = zeros(2*nnodes,1);
m(1:Ng,:) = Inpt.*meta.m;
cnt2 = 0;
cov = zeros(2*nnodes,2*nnodes);
cov2 = meta.cov;
cov2(:,Inpt==0)=0;
cov2(Inpt==0,:) = 0;
cov(1:Ng,1:Ng) = cov2;
Hi = eye(2*nnodes);

for i=1:meta.Nlvl,
    mi = m(meta.liptf{i});
    if cnt2,
        mi(mi==0) = eps;  %avoid zero times inf
        mil = log(mi);
        m(meta.levelf{i}) = exp(meta.Lf{i}*mil);
        %H2 = H;
        %H2(meta.levelf{i},meta.liptf{i}) 
        H = (m(meta.levelf{i})*ones(1,length(meta.liptf{i}))).*meta.Lf{i}./(ones(length(meta.levelf{i}),1)*mi');
        %H2(isnan(H2)) = 0;
        %cov = H2*cov*H2';
        C = H*cov(meta.liptf{i},meta.liptf{i})*H';
        C2 = H*cov(meta.liptf{i},meta.activef{i});
        cov(meta.levelf{i},meta.activef{i}) = C2;
        cov(meta.activef{i},meta.levelf{i}) = C2';
        cov(meta.levelf{i},meta.levelf{i}) = C;
    else
        m(meta.levelf{i},:) = meta.Lf{i}*mi;
        %H2 = H;
        %H2(meta.levelf{i},meta.liptf{i}) = meta.Lf{i};
        H = meta.Lf{i};  
        %cov = H2*cov*H2';
        C = H*cov(meta.liptf{i},meta.liptf{i})*H';
        C2 = H*cov(meta.liptf{i},meta.activef{i});
        cov(meta.levelf{i},meta.activef{i}) = C2;
        cov(meta.activef{i},meta.levelf{i}) = C2';
        cov(meta.levelf{i},meta.levelf{i}) = C;
    end  
    cnt2 = xor(cnt2,1);
end

m(2*nnodes) = 1;
cnt2 = mod(meta.Nlvl+1,2);
for i=meta.Nlvl:-1:1,
    mi = m(meta.liptb{i});
    if cnt2,
        mi(mi==0) = eps;  %avoid zero times inf
        mil = log(mi);
        m2 = exp(meta.Lb{i}*mil);
        m(meta.levelb{i}) = m(meta.levelb{i})+meta.Lb2{i}*m2;
        %m2 = exp(meta.Lb{i}*mil);
        H3 = meta.Lb2{i}*((m2*ones(1,length(meta.liptb{i}))).*meta.Lb{i}./(ones(size(m2))*mi'));
        H2 = Hi;
        H = H2(meta.levelb{i},meta.liptb{i}) + H3;
        %cov = H2*cov*H2';   
        C = H*cov(meta.liptb{i},meta.liptb{i})*H';
        C2 = H*cov(meta.liptb{i},meta.activeb{i});
        cov(meta.levelb{i},meta.activeb{i}) = C2;
        cov(meta.activeb{i},meta.levelb{i}) = C2';
        cov(meta.levelb{i},meta.levelb{i}) = C;
    else
        m(meta.levelb{i},:) = m(meta.levelb{i},:)+meta.Lb{i}*mi;
        H2 = Hi;
        %H2(meta.levelb{i},meta.liptb{i}) = H2(meta.levelb{i},meta.liptb{i})+meta.Lb{i};
        H = H2(meta.levelb{i},meta.liptb{i})+meta.Lb{i};
        %cov = H2*cov*H2';
        C = H*cov(meta.liptb{i},meta.liptb{i})*H';
        C2 = H*cov(meta.liptb{i},meta.activeb{i});
        cov(meta.levelb{i},meta.activeb{i}) = C2;
        cov(meta.activeb{i},meta.levelb{i}) = C2';
        cov(meta.levelb{i},meta.levelb{i}) = C;
    end
    cnt2 = xor(cnt2,1);
end

denom = m(nnodes);
vard = cov(nnodes,nnodes);
p = NaN*zeros(N,max(card));
v = NaN*zeros(N,max(card));
for i=1:N,
    for j=1:card(i),
        loc = find(meta.varmap==i);
        loc2 = loc(meta.valuemap(i,loc)==j);
        loc3 = loc2+nnodes;
        ne = length(loc2);
        m2 = zeros(ne,1);
        cov2 = zeros(ne,ne);
        cross = zeros(ne,2*nnodes);
        for k=1:ne,
            x = loc2(k);
            y = loc3(k);
            m2(k) = m(x)*m(y); %+cov(x,y);
            cov2(k,1:k-1) = m(x)*cross(1:k-1,y)'+m(y)*cross(1:k-1,x)';
            cov2(1:k-1,k) = cov2(k,1:k-1)';
            % Commented out the last two terms below to have a linear approximation for a Guassian fit
            cov2(k,k) = m(x)^2*cov(y,y)+m(y)^2*cov(x,x)+2*m(x)*m(y)*cov(x,y); %+cov(x,x)*cov(y,y)+cov(x,y)^2;
            cross(k,:) = m(x)*cov(y,:)+m(y)*cov(x,:);
        end
        mnum = sum(m2);
        var = sum(cov2(:));
        cross = sum(cross(:,nnodes));
        p(i,j) = mnum/denom;
        v(i,j) = 1/denom^2*var+mnum^2/denom^4*vard-2*mnum/denom^3*cross;
    end
end