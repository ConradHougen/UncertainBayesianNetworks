function [mo,covo] = forwardprop_so2(wcond,mm,covm)

%Calculate second order moments avoiding Taylor approx


cardo = size(wcond,2)-1;

mcond = wcond*[eye(cardo); 1/cardo*ones(1,cardo)];

np = length(mm);
card = zeros(np,1);
for i=1:np,
    card(i) = length(mm{i});
end
Ns = prod(card);
v1 = zeros(1,np);
v2 = zeros(1,np);

mo = zeros(1,cardo);
mv1 = zeros(1,np);
mv2 = zeros(1,np);
covo = zeros(cardo,cardo);


for i=1:Ns,
    for j=1:np,
        mv1(j) = mm{j}(v1(j)+1);  
    end
    val = prod(mv1);
    mo = mo+mcond(i,:)*val;
    for ii=1:Ns,
        for j=1:np,
            mv2(j) = covm{j}(v1(j)+1,v2(j)+1)+mm{j}(v2(j)+1)*mv1(j);
        end
        val = prod(mv2);
        if ii==i,
            s = cardo/wcond(i,end);
            covo = covo+(diag(mcond(i,:))+s*mcond(i,:)'*mcond(i,:))/(s+1)*val;
        else
            covo = covo+mcond(i,:)'*mcond(ii,:)*val;
        end
        for j=np:-1:1,
            v2(j) = v2(j)+1;
            if v2(j)<card(j),
                break
            else
                v2(j)=0;
            end
        end
    end
    
    for j=np:-1:1,
        v1(j) = v1(j)+1;
        if v1(j)<card(j),
            break;
        else
            v1(j) = 0;
        end
    end
end

covo = covo-mo'*mo;
    
    