function node = inferSOloopy(node, obs, NITERS)

% If niters is 0, use the MAXITERS/DELTAFLOOR stopping condition instead

MAXITERS = 100;
DELTAFLOOR = 1e-8;
N = length(node);
card = zeros(1,N);

for l=1:N
    card(l) = node(l).cardinality;
    node(l).val = [];
    node(l).childmsg = cell(length(node(l).children),1);
    node(l).childmsgc = cell(length(node(l).children),1); 
    node(l).parentmsg = cell(length(node(l).parents),1);
    node(l).parentmsgc = cell(length(node(l).parents),1);
    node(l).forward = [];
    node(l).forwardc = [];
    node(l).back = [];
    node(l).backc = [];
    node(l).delta_me = 0;
    node(l).me = zeros(1,card(l));
    
    % Init the childmsg and parentmsg fields
    % Init the childmsgc and parentmsgc fields
    nc = length(node(l).children);
    for j = 1:nc
        node(l).childmsg{j} = ones(1,card(l))/card(l);
        node(l).childmsgc{j} = zeros(card(l),card(l));
    end
    np = length(node(l).parents);
    for j = 1:np
        pcrd = card(node(l).parents(j));
        node(l).parentmsg{j} = ones(1,pcrd);
        node(l).parentmsgc{j} = zeros(pcrd,pcrd);
    end
    
end

for l=1:size(obs,1)
    node(obs(l,1)).val = obs(l,2);
end

niter = 0;
maxdelta = Inf;
while 1
    % Check stopping conditions
    if NITERS == 0
        if (maxdelta <= DELTAFLOOR || niter >= MAXITERS)
            break;
        end
    elseif NITERS > 0
        if niter >= NITERS
            break;
        end
    else
        fprintf(1, "Error: NITERS input invalid");
        return;
    end

    active = [1:N];
    id = active(randi(length(active)));
    
    % Loop through nodes in BN, sending the parent/child messages
    while ~isempty(active)
        % Send parent messages
        if isempty(node(id).val)
            if isempty(node(id).childmsg)
                crd = card(id);
                node(id).back = ones(1,crd)./crd;
                node(id).backc = zeros(crd,crd);
            else
                [node(id).back, node(id).backc] = fusion_so(node(id).childmsg,node(id).childmsgc);
                % DEBUG
%                 if id == 3
%                     fprintf(1,"node(3).back=%1.4f,$1.4f\n",node(id).back(1),node(id).back(2));
%                 end
            end
        else
            crd = card(id);
            node(id).back = zeros(1,crd);
            node(id).back(node(id).val) = 1;
            %node(id).back(node(id).val+1) = 1;
            node(id).backc = zeros(crd,crd);
        end
        for sendid = 1:length(node(id).parents)
            idx = node(id).parents(sendid);
            loc = find(node(idx).children == id);
            [node(idx).childmsg{loc},node(idx).childmsgc{loc}] = backprop_so(node(id).w,node(id).parentmsg,node(id).parentmsgc,card(node(id).parents),sendid,node(id).back,node(id).backc);
        end

        % Send child messages
        if isempty(node(id).val)
            [node(id).forward,node(id).forwardc] = forwardprop_so(node(id).w,node(id).parentmsg,node(id).parentmsgc);
            %DEBUG
%             if id == 3
%                 fprintf(1,"node(3).forward=%1.4f,%1.4f\n",node(id).forward(1),node(id).forward(2));
%             end
        else
%             m = node(id).w*[eye(card(id));1/card(id)*ones(1,card(id))];
%             node(id).forward = m;
%             s = card(id)/node(id).w(end);
%             node(id).forwardc = (diag(m)-m'*m)/(s+1);
            node(id).forward = zeros(1,card(id));
            node(id).forward(node(id).val) = 1;
            node(id).forwardc = zeros(card(id),card(id));
        end
        for sendid = 1:length(node(id).children)
            idx = node(id).children(sendid);
            loc = find(node(idx).parents==id);
            nc = length(node(id).childmsg);
            m_in = cell(nc,1);
            cov_in = cell(nc,1);
            m_in{1} = node(id).forward;
            cov_in{1} = node(id).forwardc;
            cnt = 1;
            for jj=[1:sendid-1 sendid+1:nc]
                cnt = cnt+1;
                m_in{cnt} = node(id).childmsg{jj};
                cov_in{cnt} = node(id).childmsgc{jj};
            end
            [node(idx).parentmsg{loc},node(idx).parentmsgc{loc}]= fusion_so(m_in,cov_in);
        end

        % Remove id from active list and get the next id
        loc = active==id;
        active(loc) = [];
        if ~isempty(active)
            id = active(randi(length(active)));
        end
    end
    
    % Loop through nodes and update me and mec values
    maxdelta = 0;
    for id=1:N
        if isempty(node(id).forward)
            [node(id).forward,node(id).forwardc] = forwardprop_so(node(id).w,node(id).parentmsg,node(id).parentmsgc);
        end
        if isempty(node(id).back)
            [node(id).back, node(id).backc] = fusion_so(node(id).childmsg,node(id).childmsgc);
        end
        prev_me = node(id).me;
        
        m_in = cell(2,1);
        cov_in = cell(2,1);
        m_in{1} = node(id).forward;
        cov_in{1} = node(id).forwardc;
        m_in{2} = node(id).back;
        cov_in{2} = node(id).backc;
        [node(id).me,node(id).mec] = fusion_so(m_in,cov_in);
        
        node(id).delta_me = max(abs(prev_me - node(id).me));
        
        maxdelta = max(maxdelta, node(id).delta_me);
    end
    
    niter = niter + 1;
end

%fprintf(1,"niter=%d\n",niter);

end
    
    
