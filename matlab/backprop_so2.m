function [mo,covo] = backprop_so2(wcond,mm,covm,card,parent,ml,covl)

np = length(mm);
% card = zeros(np+1,1);
% for i=1:np,
%     card(i) = length(mm{i});
% end
% card(end) = length(ml);

%card = [card length(ml)];

cardl = length(ml);

mcond = wcond*[eye(cardl); 1/card(end)*ones(1,cardl)];

%relparents = setdiff(1:np+1,parent);
relparents = setdiff(1:np,parent);

Ns = prod(card);


%np1 = np+1;
m = zeros(1,np);
m2 = zeros(1,np);

mo = zeros(1,card(parent));

cov = zeros(card(parent),card(parent));



vv = zeros(1,np);
vv2 = zeros(1,np);

for i=1:Ns,
    v = vv+1;
    %k1 = floor((i-1)/card(end))+1;
    for j=[1:parent-1 parent+1:np],
        m(j) = mm{j}(v(j));
    end
    %m(end) = ml(v(end));
    
    mp = prod(m(relparents));
    mo(v(parent)) = mo(v(parent))+sum(ml.*mcond(i,:)*mp);
    
    for ii=1:Ns,
        v2 = vv2+1;
        %k2 = floor((ii-1)/card(end))+1;
        for j=[1:parent-1 parent+1:np],
            m2(j) = covm{j}(v(j),v2(j))+mm{j}(v(j))*mm{j}(v2(j));
        end
        %m2(end) = covl(v(end),v2(end))+ml(v(end))*ml(v2(end));
        mp = prod(m2(relparents));
        if (ii==i),
            s = card(end)/wcond(i,end);
            cov(v(parent),v2(parent)) = cov(v(parent),v2(parent))+sum(sum((covl+ml'*ml).*(diag(mcond(i,:))+s*mcond(i,:)'*mcond(ii,:))))/(s+1)*mp;
        else
            cov(v(parent),v2(parent)) = cov(v(parent),v2(parent)) + sum(sum((covl+ml'*ml).*(mcond(i,:)'*mcond(ii,:))))*mp;
        end
        for j=np:-1:1,
            vv2(j) = vv2(j)+1;
            if vv2(j)==card(j),
                vv2(j) = 0;
            else
                break
            end
        end   
    end

    for j=np:-1:1,
        vv(j) = vv(j)+1;
        if vv(j)==card(j),
            vv(j) = 0;
        else
            break
        end
    end           
end

cov = cov - mo'*mo;

denom = sum(mo);
mo = mo/denom;

H = (eye(card(parent))-mo'*ones(1,card(parent)))/denom;

covo = H*cov*H';