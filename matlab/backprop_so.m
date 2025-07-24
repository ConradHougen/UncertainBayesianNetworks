function [mo,covo] = backprop_so(wcond,mm,covm,card,parent,ml,covl)

np = length(mm);
% card = zeros(np+1,1);
% for i=1:np,
%     card(i) = length(mm{i});
% end
% card(end) = length(ml);

card = [card length(ml)];

mcond = wcond*[eye(card(end)); 1/card(end)*ones(1,card(end))];

relparents = setdiff(1:np+1,parent);

H = cell(np,1);

for i=1:np,
    H{i} = zeros(card(parent),card(relparents(i)));
end
Ns = prod(card);

H2 = zeros(card(parent),Ns);


np1 = np+1;
m = zeros(1,np1);

mo = zeros(1,card(parent));

Ns2 = Ns/card(end);

covmc = zeros(Ns2,Ns2);

vv = zeros(1,np1);

for i=1:Ns,
    v = vv+1;
    if vv(end)==0,
        ii = (i-1)/card(end);
        s = card(end)/wcond(ii+1,end);
        covmc(ii*card(end)+1:(ii+1)*card(end),ii*card(end)+1:(ii+1)*card(end)) = (diag(mcond(ii+1,:))-mcond(ii+1,:)'*mcond(ii+1,:))/(s+1);
        ii = ii+1;
    end
    for j=[1:parent-1 parent+1:np],
        m(j) = mm{j}(v(j));
    end
    m(end) = ml(v(end));
    
    mp = [1 cumprod(m(relparents(1:end-1)))];
  
    for j=np:-1:1,
        mf = mp(j)*prod(m(relparents(j+1:np)));
        jj = relparents(j); 
        %idp = setdiff(relparents,jj);
        %H{j}(v(parent),v(jj)) = H{j}(v(parent),v(jj)) + prod(m(idp))*mcond(ii,v(end));
        H{j}(v(parent),v(jj)) = H{j}(v(parent),v(jj)) + mf*mcond(ii,v(end));
    end
    mp = mf*m(relparents(1)); %prod(m(relparents));
    H2(v(parent),i) = mp;
    mo(v(parent)) = mo(v(parent))+mcond(ii,v(end))*mp;
    for j=np1:-1:1,
        vv(j) = vv(j)+1;
        if vv(j)==card(j),
            vv(j) = 0;
        else
            break
        end
    end           
end

cov = zeros(card(parent),card(parent));

for i=1:np-1,
    j = relparents(i);
    cov = cov+H{i}*covm{j}*H{i}';
end

cov = cov+H{np}*covl*H{np}'+H2*covmc*H2';

denom = sum(mo);
mo = mo/denom;

H = (eye(card(parent))-mo'*ones(1,card(parent)))/denom;

covo = H*cov*H';