% cdh_net2circuit
% Translate Bayes Net into sum-product network
function [circ, meta] = cdh_net2circuit(net)

n = length(net);
card = zeros(1,n);
active = zeros(1,n);

for i=1:n
    card(i) = net(i).cardinality;
    active(i) = isempty(net(i).children); % set active nodes as leaf nodes
end

active = find(active == 1); % Translate into indices of activ enodes
nground = 0; % number of ground (leaf) nodes in sum product network
sz = zeros(n,1);

% accumulate
for i=1:n
    % sz(i) is the number of cond probabilities associated with node i
    sz(i) = prod([card(i) card(net(i).parents)]);
    nground = nground + sz(i);
end

p = zeros(nground, 1); % conditional probabilities
varmap = zeros(1, nground);
valuemap = zeros(n, nground);

cnt = 0;
circ = struct('fun',cell(nground,1),'parents',cell(nground,1),'children',cell(nground,1));

for i=1:n
    idx = net(i).parents;
    np = length(idx);
    val = zeros(1,np);
    cd = card(idx);

    % If node has no parents
    if np==0
        for jj=1:card(i)
            cnt = cnt+1;
            circ(cnt).fun = 'ground';
            varmap(1,cnt) = i;
            valuemap(i,cnt) = jj;
            p(cnt) = net(i).p(jj,1);

        end

    % Node i has parents
    else
        sz2 = sz(i)/card(i);
        for j=1:sz2
            for jj=1:card(i)
                cnt = cnt+1;
                circ(cnt).fun = 'ground';
                varmap(1,cnt) = i;
                valuemap(idx,cnt) = val'+1;
                valuemap(i,cnt) = jj;
                p(cnt) = net(i).p(jj,j);
            end

            val(np) = val(np)+1;
            for k=np:-1:2
                if val(k) == cd(k)
                    val(k) = 0;
                    val(k-1) = val(k-1)+1;
                else
                    break
                end
            end
        end
    end
end

opennodes = zeros(1,0);

while ~isempty(active)
    i = active(1);
    loc = find(varmap==i);
    idx = net(i).parents;
    N = length(idx);
    val = zeros(1,N);
    cd = card(idx);

    locm = find((valuemap(i,nground+1:end)>0)&opennodes)+nground;
    idx = setdiff(find(sum(valuemap(:,[loc locm]),2)>1),i);

    if ~isempty(locm)
        N = length(idx);
        val = zeros(1,N);
        cd = card(idx);
        sz2 = prod(cd);
        loc = union(loc,locm); 
        sz3 = length(loc);
        for j=1:sz2                    
            loc2 = loc(sum((valuemap(idx,loc)==0)+(valuemap(idx,loc)==(val+1)'*ones(1,sz3)),1)==N);

            for ii=1:card(i)
                loc3 = find(valuemap(i,loc2)==ii);
                cnt = cnt+1;
                circ(cnt).fun = 'and';
                circ(cnt).children = loc2(loc3);
                opennodes(cnt-nground) = 0;
                valuemap([i; idx],cnt) = [ii; (val+1)'];
                for jj=loc2(loc3),
                    circ(jj).parents = [circ(jj).parents cnt];
                    if jj>nground,
                        opennodes(jj-nground) = 0;
                    end
                end
            end
            for jj=cnt-card(i)+1:cnt,
                circ(jj).parents = cnt+1;
            end
            cnt = cnt+1;
            circ(cnt).fun = 'or';
            circ(cnt).children = cnt-card(i):cnt-1;
            valuemap(idx,cnt) = (val+1)';
            opennodes(cnt-nground) = 1;
            if sz2==1,
                break,
            end
            val(N) = val(N)+1;
            for k=N:-1:2,
                if val(k) == cd(k),
                    val(k) = 0;
                    val(k-1) = val(k-1)+1;
                else
                    break
                end
            end    
        end
    elseif N==0,
        cnt = cnt+1;
        circ(cnt).fun = 'or';
        circ(cnt).children = loc;
        valuemap(:,cnt) = zeros(n,1);
        for j=loc,
            circ(loc).parents = [circ(loc).parents cnt];
        end
        opennodes(cnt-nground) = 0;
    else
        sz2 = sz(i)/card(i);
        for j=1:sz2,
            cnt = cnt+1;
            circ(cnt).fun = 'or';
            loc2 = find(sum(valuemap(idx,loc)==(val+1)'*ones(1,sz(i)),1)==N);
            circ(cnt).children = loc(loc2);
            valuemap(idx,cnt) = (val+1)';
            opennodes(cnt-nground) = 1;
            for jj=loc2,
                circ(loc(jj)).parents = [circ(loc(jj)).parents cnt];
            end
            val(N) = val(N)+1;
            for k=N:-1:2,
                if val(k) == cd(k),
                    val(k) = 0;
                    val(k-1) = val(k-1)+1;
                else
                    break
                end
            end
        end
    end
    active = active(2:end);
    for ii=net(i).parents,
        net(ii).children = setdiff(net(ii).children,i);
        if isempty(net(ii).children),
            active = union(active,ii);
        end
    end
end

meta.varmap = varmap;
meta.valuemap = valuemap;
meta.p = p;
meta.opennodes = opennodes;
    
end