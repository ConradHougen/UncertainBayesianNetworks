function [p,m,active,circ] = cdh_inferProb(circ,meta,obs)

m = meta.p;
varmap = meta.varmap;
valuemap = meta.valuemap;

nground = length(m);

nnodes = length(circ);
m = [m;zeros(nnodes-nground,1)];

for i=1:nnodes
    circ(i).receive = zeros(size(circ(i).children));
end

active = [];
for i=1:nground
    for j=circ(i).parents
        circ(j).receive(circ(j).children==i)=1;
        if isempty(circ(j).receive(circ(j).receive==0))
            active = union(active,j);
        end
    end
end

nobs = size(obs,1);
indicators = ones(nground, 1);
for i=1:nobs
    loc = find(varmap==obs(i,1));
    loc = loc(valuemap(obs(i,1),loc)~=obs(i,2));
    m(loc) = 0;
    indicators(loc) = 0;
end

% Forward pass
while ~isempty(active)
    i = active(1);
    loc = circ(i).children;
    if circ(i).fun(1) == 'o'
        m(i) = sum(m(loc));
    elseif circ(i).fun(1) == 'a'
        m(i) = prod(m(loc));
    end
    for j=circ(i).parents
        circ(j).receive(circ(j).children==i)=1;
        if isempty(circ(j).receive(circ(j).receive==0))
            active = union(active,j);
        end
    end
    active = setdiff(active,i);
end

denom = m(end);

active = nnodes;

m = [m; zeros(nnodes,1)];
m(end) = 1;

for i=1:nnodes
    circ(i).receive = zeros(size(circ(i).parents));
end

% Backward pass
while ~isempty(active)
    i = active(1);
    if circ(i).fun(1) == 'o'      
        for j=circ(i).children
            m(j+nnodes) = m(j+nnodes)+m(i+nnodes);
            circ(j).receive(circ(j).parents==i)=1;
            if isempty(circ(j).receive(circ(j).receive==0))
                active = union(active,j);
            end
        end   
    elseif circ(i).fun(1) == 'a'
        pval = [1; cumprod(m(circ(i).children(1:end-1)))];
        for ii=length(circ(i).children):-1:1
            j=circ(i).children(ii);
            sendmsg =prod([pval(ii); m(circ(i).children(ii+1:end))])*m(i+nnodes);
            m(j+nnodes) = m(j+nnodes) + sendmsg;
            circ(j).receive(circ(j).parents==i)=1;
            if isempty(circ(j).receive(circ(j).receive==0))
                active = union(active,j);
            end
        end
    end
    active = setdiff(active,i);
end

% Set partial derivatives to zero as per indicator nodes
%m(nnodes+1:nnodes+nground) = m(nnodes+1:nnodes+nground) .* (m(1:nground) > 0);
m(nnodes+1:nnodes+nground) = m(nnodes+1:nnodes+nground) .* indicators;


n = size(valuemap,1);
card = max(valuemap,[],2);
p = NaN*zeros(n,max(card));
for i=1:n
    for j=1:card(i)
        loc = find(varmap==i);
        loc2 = loc(valuemap(i,loc)==j);
        p(i,j) = sum(m(loc2).*m((loc2)+nnodes))/denom; %p(i,j) = p(node i = value j | evidence)
    end
end

end