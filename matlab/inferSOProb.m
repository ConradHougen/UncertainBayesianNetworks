function [p,v,m,cov,active,circ] = inferSOProb(circ,meta,obs)

m = meta.m;
cov = meta.cov;
varmap = meta.varmap;
valuemap = meta.valuemap;

nground = length(m);

nnodes = length(circ);

sz = 2*nnodes-nground;
m = [m;zeros(sz,1)];
cov = [cov zeros(nground,sz); zeros(sz,nground) zeros(sz,sz)];

for i=1:nnodes,
    circ(i).receive = zeros(size(circ(i).children));
end

active = [];
for i=1:nground,
    for j=circ(i).parents,
        circ(j).receive(circ(j).children==i)=1;
        if isempty(circ(j).receive(circ(j).receive==0)),
            active = union(active,j);
        end
    end
end

nobs = size(obs,1);
for i=1:nobs,
    loc = find(varmap==obs(i,1));
    loc = loc(valuemap(obs(i,1),loc)~=obs(i,2));
    m(loc) = 0;
    cov(loc,:) = 0;
    cov(:,loc) = 0;
end

while ~isempty(active),
    i = active(1);
    loc = circ(i).children;
    if circ(i).fun(1) == 'o',
        m(i) = sum(m(loc));
        cov(i,:) = sum(cov(loc,:),1);
        cov(:,i) = cov(i,:)';
        cov(i,i) = sum(sum(cov(loc,loc)));
    elseif circ(i).fun(1) == 'a',
        mf = prod(m(loc));
        if mf==0,
            m(i) = 0;
            cov(i,:) = 0;
            cov(:,i) = 0;
        else
            o = mf./m(loc);
            c2 = o'*cov(loc,loc)*o;
            c = cov(:,loc)*o;
            cov(:,i) = c;
            cov(i,:) = c';
            cov(i,i) = c2;
            m(i) = mf;
        end
%         idx = loc(1);
%         m(i) = m(idx);
%         cov(:,i) = cov(:,idx);
%         cov(i,:) = cov(idx,:);
%         cov(i,i) = cov(idx,idx);
%         for idx = loc(2:end),
%             c = m(idx)*cov(:,i)+m(i)*cov(:,idx);
%             % Commented out the last two terms below to have a linear approximation for a Guassian fit
%             c2 = m(idx)^2*cov(i,i)+m(i)^2*cov(idx,idx)+2*m(idx)*m(i)*cov(idx,i); %+cov(idx,idx)*cov(i,i)+cov(idx,i)^2;
%             m(i) = m(i)*m(idx); %+cov(idx,i);
%             cov(:,i) = c;
%             cov(i,:) = c';
%             cov(i,i) = c2;
%         end
        %%%%%%%%%%%%%%%%%%%%%Start here%%%%%%%%%%%%%%%%%%%%%%% 
%         m(i) = prod(m(loc));
    end
    for j=circ(i).parents,
         circ(j).receive(circ(j).children==i)=1;
        if isempty(circ(j).receive(circ(j).receive==0)),
            active = union(active,j);
        end
    end
    active = setdiff(active,i);
end

denom = m(nnodes);
vard = cov(nnodes,nnodes);

active = nnodes;

m(end) = 1;

for i=1:nnodes,
    circ(i).receive = zeros(size(circ(i).parents));
end


while ~isempty(active),
    i = active(1);
    if circ(i).fun(1) == 'o',        
        for j=circ(i).children,     
            m(j+nnodes) = m(j+nnodes)+m(i+nnodes);
            c = cov(:,j+nnodes)+cov(:,i+nnodes);
            c2 = c(j+nnodes)+c(i+nnodes);
            cov(:,j+nnodes) = c;
            cov(j+nnodes,:) = c';
            cov(j+nnodes,j+nnodes) = c2;
            circ(j).receive(circ(j).parents==i)=1;
            if isempty(circ(j).receive(circ(j).receive==0)),
                active = union(active,j);
            end
        end   
    elseif circ(i).fun(1) == 'a'
        nc = length(circ(i).children);
        mf = m(i+nnodes);
        mc = m(circ(i).children);
        lz = find(mc==0);
        if mf>0,
            lng = length(lz);
            if lng==0,
                mp = prod(mc);
                mf = mf*mp;
                loc = [circ(i).children i+nnodes circ(i).children+nnodes];
                o = [(mf./mc)*(1./mc)'-diag(mf./mc.^2) mp./mc eye(nc)];
                c2 = o*cov(loc,loc)*o';
                c = cov(:,loc)*o';
                loc = circ(i).children+nnodes;
                m(loc) = m(loc)+mf./m(loc-nnodes);
                cov(:,loc) = c;
                cov(loc,:) = c';
                cov(loc,loc) = c2;     
            elseif lng==1, 
                mc = mc([1:lz-1 lz+1:end]);
                loc = [circ(i).children([1:lz-1 lz+1:end]) i+nnodes circ(i).children(lz)+nnodes];
                mp = prod(mc);
                mf = mf*mp;
                o = [mf./mc' mp 1];
                c2 = o*cov(loc,loc)*o';
                c = cov(:,loc)*o';
                loc = circ(i).children(lz)+nnodes;
                m(loc) = m(loc)+mf;
                cov(:,loc) = c;
                cov(loc,:) = c';
                cov(loc,loc) = c2;
            end
        end     
        
%         m2 = zeros(nc,1);
%         cov2 = zeros(nc,nc);
%         cross = zeros(nc,2*nnodes);
%         %idx = circ(i).children(1);
%         m2(1) = 1; %m(idx);
%         %cov2(1,1) = cov(idx,idx);
%         %cross(1,:) = cov(idx,:); %[circ(i).children i+nnodes]);
%         for j=2:nc,
%             idx = circ(i).children(j-1);
%             m2(j) = m2(j-1)*m(idx); %+cross(j-1,idx);
%             if m2(j)==0, %avoid unnecessary computations when the product is zero
%                  break;
%             end
%             cross(j,:) = m2(j-1)*cov(idx,:)+m(idx)*cross(j-1,:);
%             cov2(j,:) = m2(j-1)*cross(:,idx)'+m(idx)*cov2(j-1,:);
%             cov2(:,j) = cov2(j,:)';
%             % Commented out the last two terms below to have a linear approximation for a Guassian fit
%             cov2(j,j) = m2(j-1)^2*cov(idx,idx)+m(idx)^2*cov2(j-1,j-1)+2*m2(j-1)*m(idx)*cross(j-1,idx); %+cov(idx,idx)*cov2(j-1,j-1)+cross(j-1,idx)^2;
%         end
%         for j=nc:-1:1,
%             for k=j+1:nc,
%                 idx = circ(i).children(k);
%                 m2(k) = m2(k-1)*m(idx); %+cross(k-1,idx);
%                 if m2(k)==0,  %avoid unnecessary computations when the product is zero
%                     cross(k:nc,:) = 0;
%                     cov2(k:nc,:) = 0;
%                     cov2(:,k:nc) = 0;
%                     break;
%                 end
%                 cross(k,:) = m2(k-1)*cov(idx,:)+m(idx)*cross(k-1,:);
%                 cov2(k,:) = m2(k-1)*cross(:,idx)'+m(idx)*cov2(k-1,:);
%                 cov2(:,k) = cov2(k,:)';
%                 % Commented out the last two terms below to have a linear approximation for a Guassian fit
%                 cov2(k,k) = m2(k-1)^2*cov(idx,idx)+m(idx)^2*cov2(k-1,k-1)+2*m2(k-1)*m(idx)*cross(k-1,idx); %+cov(idx,idx)*cov2(k-1,k-1)+cross(k-1,idx)^2;
%             end            
%             msend = m2(end)*m(i+nnodes); %+cross(end,i+nnodes);
%             crsend = m2(end)*cov(i+nnodes,:)+m(i+nnodes)*cross(end,:);
%             % Commented out the last two terms below to have a linear approximation for a Guassian fit
%             csend = m2(end)^2*cov(i+nnodes,i+nnodes)+m(nnodes+i)^2*cov2(end,end)+2*m2(end)*m(nnodes+i)*cross(end,nnodes+i); %+cov2(end,end)*cov(i+nnodes,i+nnodes)+cross(end,nnodes+i)^2;
%             idx = circ(i).children(j)+nnodes;
%             m(idx) = m(idx) + msend;         
%             var = cov(idx,idx);
%             cov(idx,:) = cov(idx,:) + crsend;
%             cov(:,idx) = cov(idx,:)';         
%             cov(idx,idx) = var+csend+2*crsend(idx);
        for j=1:nc,
            idx = circ(i).children(j);
            circ(idx).receive(circ(idx).parents==i)=1;
            if isempty(circ(idx).receive(circ(idx).receive==0)),
                active = union(active,idx);
            end
        end
    end
    active = setdiff(active,i);
end

n = size(valuemap,1);
card = max(valuemap,[],2);
p = NaN*zeros(n,max(card));
v = NaN*zeros(n,max(card));
for i=1:n,
    for j=1:card(i),
        loc = find(varmap==i);
        loc2 = loc(valuemap(i,loc)==j);
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









