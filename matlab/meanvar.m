function [pc,var] = meanvar(circ,meta,obs)

m = meta.m;
varmap = meta.varmap;
valuemap = meta.valuemap;

n = size(valuemap,1);
card = max(valuemap,[],2);
maxcard = max(card);
pc = NaN*ones(n,maxcard);
var = NaN*ones(n,maxcard);

nground = length(m);

diff = zeros(nground,maxcard-1);

nnodes = length(circ);
m = [m;zeros(nnodes-nground,1)];

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
zloc = [];
for i=1:nobs,
    loc = find(varmap==obs(i,1));
    loc = loc(valuemap(obs(i,1),loc)~=obs(i,2));
    m(loc) = 0;
    pc(obs(i,1),1:card(obs(i,1))) = 0;
    var(obs(i,1),1:card(obs(i,1))) = 0;
    pc(obs(i,1),obs(i,2)) = 1;
    zloc = union(zloc,loc);
end

query = setdiff(1:n,obs(:,1));


while ~isempty(active),
    i = active(1);
    loc = circ(i).children;
    if circ(i).fun(1) == 'o', 
        m(i) = sum(m(loc));
    elseif circ(i).fun(1) == 'a',
        m(i) = prod(m(loc));
    end
    for j=circ(i).parents,
         circ(j).receive(circ(j).children==i)=1;
        if isempty(circ(j).receive(circ(j).receive==0)),
            active = union(active,j);
        end
    end
    active = setdiff(active,i);
end

denom = m(end);

active = nnodes;

m = [m; zeros(nnodes,1)];
m(end) = 1;

for i=1:nnodes,
    circ(i).receive = zeros(size(circ(i).parents));
end


while ~isempty(active),
    i = active(1);
    if circ(i).fun(1) == 'o',        
        for j=circ(i).children,
            m(j+nnodes) = m(j+nnodes)+m(i+nnodes);
            circ(j).receive(circ(j).parents==i)=1;
            if isempty(circ(j).receive(circ(j).receive==0)),
                active = union(active,j);
            end
        end   
    elseif circ(i).fun(1) == 'a',
        for j=circ(i).children,
            sendmsg =prod(m(setdiff(circ(i).children,j)))*m(i+nnodes);
            m(j+nnodes) = m(j+nnodes) + sendmsg;
            circ(j).receive(circ(j).parents==i)=1;
            if isempty(circ(j).receive(circ(j).receive==0)),
                active = union(active,j);
            end
        end
    end
    active = setdiff(active,i);
end

msave = m;

dpe = m(nnodes+1:nnodes+nground);
dpe(zloc) = 0;


for q = query,
    for qv = 1:card(q)-1,
        
        m = msave;
        
        loc = find(varmap==q);
        loc = loc(valuemap(q,loc)~=qv);
        m(loc) = 0;

        zloc2 = union(zloc,loc);

        m(nground+1:end) = 0;

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

        while ~isempty(active),
            i = active(1);
            loc = circ(i).children;
            if circ(i).fun(1) == 'o',
                m(i) = sum(m(loc));
            elseif circ(i).fun(1) == 'a',
                m(i) = prod(m(loc));
            end
            for j=circ(i).parents,
                circ(j).receive(circ(j).children==i)=1;
                if isempty(circ(j).receive(circ(j).receive==0)),
                    active = union(active,j);
                end
            end
            active = setdiff(active,i);
        end

        num = m(nnodes);

        active = nnodes;

        %m = [m; zeros(nnodes,1)];
        m(end) = 1;

        for i=1:nnodes,
            circ(i).receive = zeros(size(circ(i).parents));
        end


        while ~isempty(active),
            i = active(1);
            if circ(i).fun(1) == 'o',        
                for j=circ(i).children,
                    m(j+nnodes) = m(j+nnodes)+m(i+nnodes);
                    circ(j).receive(circ(j).parents==i)=1;
                    if isempty(circ(j).receive(circ(j).receive==0)),
                        active = union(active,j);
                    end
                end   
            elseif circ(i).fun(1) == 'a',
                for j=circ(i).children,
                    sendmsg =prod(m(setdiff(circ(i).children,j)))*m(i+nnodes);
                    m(j+nnodes) = m(j+nnodes) + sendmsg;
                    circ(j).receive(circ(j).parents==i)=1;
                    if isempty(circ(j).receive(circ(j).receive==0)),
                        active = union(active,j);
                    end
                end
            end
            active = setdiff(active,i);
        end

        pc(q,qv) = num/denom;



        dphe = m(nnodes+1:nnodes+nground);
        dphe(zloc2) = 0;

        dtheta = 1/denom*(dphe-pc(q,qv)*dpe);
        
        diff(:,qv) = dtheta;

        % n = size(valuemap);
        % card = max(valuemap,[],2);
        % var = 0;
        % 
        % 
        % mout = m;
        % c = diag(meta.cov);
        % m = meta.m;
        % for i=1:n,
        %    loc = find(varmap==i);
        %    parents = setdiff(find(sum(valuemap(:,loc),2)>0),i);
        %    if isempty(parents),
        %        idx = loc(1);
        %        iSp1 = c(idx)/m(idx)./(1-m(idx));
        %        var = var + iSp1*(sum(dtheta(loc).^2.*m(loc))-sum(dtheta(loc).*m(loc)).^2);
        %    else
        %        nl = length(loc);
        %        np = length(parents);
        %        v = zeros(1,np);
        %        j = 1;
        %        while j>0,
        %            loc2 = loc(sum((valuemap(parents,loc)-(v+1)'*ones(1,nl))==0,1)==np);
        %            idx = loc2(1);
        %            iSp1 = c(idx)/m(idx)./(1-m(idx));
        %            var = var + iSp1*(sum(dtheta(loc2).^2.*m(loc2))-sum(dtheta(loc2).*m(loc2)).^2);
        %            j=np;
        %            while j>0,
        %                v(j) = v(j)+1;
        %                if v(j)<card(parents(j)),
        %                    break;
        %                else
        %                    v(j) = 0;
        %                    j = j-1;
        %                end
        %            end
        %        end
        %    end
        % end       
        cov = meta.cov(1:nground,1:nground);
        var(q,qv) = dtheta'*cov*dtheta;       
    end
    
    %compute variance of the final variable value from pervious variable
    %values to avoid calculations of derivatives unnecessarily
    qv = card(q);
    var(q,qv) = 0;
    pc(q,qv) = 1-sum(pc(q,1:qv-1));
    for qv1 = 1:card(q)-1,
        for qv2 = 1:card(q)-1,
            var(q,qv) = var(q,qv)+diff(:,qv1)'*cov*diff(:,qv2);
        end
    end            
end






