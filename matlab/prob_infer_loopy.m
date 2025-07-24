% Loopy BP
function node = prob_infer_loopy(node,obs)

MAXITERS = 100;
DELTAFLOOR = 0.00000001;
N = length(node);
card = zeros(N,1);

for i=1:N
    card(i) = node(i).cardinality;
    node(i).val = [];
    node(i).childmsg = ones(length(node(i).children),card(i))./card(i);
    np = max(card(node(i).parents));
    node(i).parentmsg = ones(length(node(i).parents),np);
    node(i).forward = [];
    node(i).back = [];
    node(i).pe = zeros(1,card(i));
    node(i).delta_pe = 0;
end

for i=1:size(obs,1)
    node(obs(i,1)).val = obs(i,2);    
end

% Use Belief Update algorithm 
niter = 0;
maxdelta = Inf;
while maxdelta > DELTAFLOOR && niter < MAXITERS
    %fprintf(1, "niter=%d\n", niter);
    active = [1:N];
    id = active(randi(length(active)));
    
    % Loop through all nodes in the BN, sending all parent/child messages
    while(~isempty(active))
        
        % Send parent messages
        if isempty(node(id).val)
            node(id).back = prob_fuse(node(id).childmsg);
        else
            node(id).back = zeros(1,card(id));
            node(id).back(node(id).val) = 1;
        end
        for sendid = 1:length(node(id).parents)
            idx = node(id).parents(sendid);
            loc = find(node(idx).children == id);
            node(idx).childmsg(loc,:) = prob_backward(node(id).p,node(id).parentmsg,card(node(id).parents),sendid,node(id).back);
        end
        
        % Send child messages
        if isempty(node(id).val)
            node(id).forward = prob_forward(node(id).p,node(id).parentmsg,card(node(id).parents));
        else
            node(id).forward = zeros(1,card(id));
            node(id).forward(node(id).val) = 1;
        end
        for sendid = 1:length(node(id).children)
            idx = node(id).children(sendid);
            loc = find(node(idx).parents==id);
            node(idx).parentmsg(loc,1:card(id)) = prob_fuse([node(id).forward; node(id).childmsg([1:sendid-1 sendid+1:end],:)]);
        end
        
        % Remove id node from active list and get the next id
        loc = active==id;
        active(loc) = [];
        if ~isempty(active)
            id = active(randi(length(active)));
        end
    end
    
    % Loop through all nodes and update pe values
    maxdelta = 0;
    for id=1:N
        if isempty(node(id).forward)
            node(id).forward = prob_forward(node(id).p,node(id).parentmsg,card(node(id).parents));
        end
        if isempty(node(id).back)
            node(id).back = prob_fuse(node(id).childmsg);
        end
        prev_pe = node(id).pe;
        node(id).pe = prob_fuse([node(id).forward; node(id).back]);
        node(id).delta_pe = max(abs(prev_pe - node(id).pe));
        %fprintf(1, "id=%d, pe(1)=%f, prev_pe(1)=%f\n", id, node(id).pe(1), prev_pe(1));
        
        % Update maximum observed difference in estimated pe (stopping cond.)
        maxdelta = max(maxdelta, node(id).delta_pe);
    end
            
    % Increment number of runs to hard cap to prevent infinite iterations
    niter = niter + 1;
end

%fprintf(1,"niter=%d\n",niter);

end