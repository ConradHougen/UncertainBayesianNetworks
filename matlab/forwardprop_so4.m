function [mo,covo] = forwardprop_so(wcond,mm,covm)


cardo = size(wcond,2)-1;

mcond = wcond*[eye(cardo); 1/cardo*ones(1,cardo)];

np = length(mm);
card = zeros(np,1);
for i=1:np,
    card(i) = length(mm{i});
end
Ns = prod(card);
v = zeros(1,np);

mo = zeros(1,cardo);
mv = zeros(1,np);
covo = zeros(cardo,cardo);
H = cell(np,1);
for j=1:np,
    H{j} = zeros(cardo,card(j));
end
for i=1:Ns,
    for j=1:np,
        mv(j) = mm{j}(v(j)+1);  
    end
    %mp = prod(mv);
    mp = [1 cumprod(mv(1:np-1))];
    for j=np:-1:1,
        mf = mp(j)*prod(mv(j+1:np));
%         if mv(j)==0,
%             mf = prod(mv([1:j-1 j+1:np]));
%         else
%             mf = mp/mv(j);
%         end
        H{j}(:,v(j)+1) = H{j}(:,v(j)+1)+mcond(i,:)'*mf; %prod(mv([1:j-1 j+1:np]));
    end
    val = prod(mv);
    s = cardo/wcond(i,end);
    covo = covo + val^2*(diag(mcond(i,:))-mcond(i,:)'*mcond(i,:))/(s+1);
    mo = mo + mcond(i,:)*val;
    for j=np:-1:1,
        v(j) = v(j)+1;
        if v(j)<card(j),
            break;
        else
            v(j) = 0;
        end
    end
end

for j=1:np,
    covo = covo + H{j}*covm{j}*H{j}';
end
    
    