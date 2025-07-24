function node = inferSO2(node)

N = length(node);

% active = [];

card = zeros(1,N);
for i=1:N,
    card(i) = node(i).cardinality;
end
for i=1:N,
    node(i).sendchildren = zeros(size(node(i).children));
    node(i).recchildren = zeros(size(node(i).children));
    node(i).sendparents = zeros(size(node(i).parents));
    node(i).recparents =  zeros(size(node(i).parents));
    node(i).childmsg = cell(length(node(i).children),1);
    node(i).childmsgc = cell(length(node(i).children),1); %zeros(length(node(i).children),card(i)+1);
    node(i).parentmsg = cell(length(node(i).parents),1); %zeros(length(node(i).parents),max(card(node(i).parents))+1);
    node(i).parentmsgc = cell(length(node(i).parents),1);
    node(i).forward = [];
    node(i).forwardc = [];
    node(i).back = [];
    node(i).backc = [];
%     if (length(node(i).recparents)+length(node(i).recchildren)==1),
%         active = [active i];
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Following code enables uncertain measurement at a single parent leaf node that reprents the output of
% an uncertian machine learning. that really represents the likelihood of a particular measurement

for i=1:N,
    l = length(node(i).value)-1;
    if l>0,
        if abs(sum(node(i).value)-1)>eps,
            error('Uncertain measurement is a subjective opinion that must sum to one!');
        end
        np = length(node(i).parents);
        nc = length(node(i).children);
        if (nc>0)|(np~=1),
            error(sprintf('Node %d cannot have an uncertain measurement! It is not a single parent leaf node.',i));
        end
        lp = size(node(node(i).parents).w,2)-1;
        idx = node(i).parents;
        if lp~=l,
            error(sprintf('Uncertaint measurement at Node %d must match the cardinality for the parent Node %d!',i,idx));
        end
        node(i).sendparents = 1;
        node(i).recparents = 1;
        l = size(node(i).w,2)-1;
        node(i).back = [1/l*ones(1,l)];
        node(i).backc = zeros(l,l);
        node(i).forward = node(i).back;
        node(i).forwardc = node(i).backc;
        loc = find(node(idx).children==i);
        node(idx).recchildren(loc) = 1;
        node(idx).sendchildren(loc) = 1;
        m = node(i).value*[eye(lp); 1/lp*ones(1,lp)];
        node(idx).childmsg{loc} = m;
        s = lp/node(i).value(end);
        node(idx).chilmsgm{loc} = (diag(m)-m'*m)/(s+1);
    end
end

active = [];

for i=1:N,
    if (sum(node(i).recparents==0)+sum(node(i).recchildren==0)==1),
        active = [active i];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while ~isempty(active),
    id = min(active);
    recparent = find(node(id).recparents==0);
    recchild = find(node(id).recchildren==0);
    if length(recparent)+length(recchild)==1,
        if isempty(recparent),
            sendid = recchild;
            schild = 1;
        else
            sendid = recparent;
            schild = 0;
        end
    else
        sendid = find(node(id).sendparents==0,1);
        schild = 0;
        if isempty(sendid),
            sendid = find(node(id).sendchildren==0,1);
            schild = 1;
        end
    end
    if schild==0,
        if isempty(node(id).back),
            if isempty(node(id).value),
                if isempty(node(id).childmsg),
                    l = card(id);
                    node(id).back = 1/l*ones(1,l); % 0];
                    node(id).backc = zeros(l,l);
                else
                    [node(id).back, node(id).backc] = fusion_so(node(id).childmsg,node(id).childmsgc);
                end
            else
                node(id).back = zeros(1,card(id)); %zeros(1,card(id)+1);
                node(id).back(node(id).value+1) = 1;
                %node(id).back(node(id).value+1) = 1;
                node(id).backc = zeros(card(id),card(id));
            end
        end
        
        node(id).sendparents(sendid) = 1;
        idx = node(id).parents(sendid);
        loc = find(node(idx).children==id);
        %node(idx).childmsg(loc,:) = approx_backprop(node(id).w,node(id).parentmsg,card(node(id).parents),sendid,node(id).back);

        
        [node(idx).childmsg{loc},node(idx).childmsgc{loc}] = backprop_so2(node(id).w,node(id).parentmsg,node(id).parentmsgc,card(node(id).parents),sendid,node(id).back,node(id).backc);
        
        node(idx).recchildren(loc) = 1;
        leftsend = sum(node(idx).sendparents==0)+sum(node(idx).sendchildren==0);
        leftrec = sum(node(idx).recchildren==0)+sum(node(idx).recparents==0);
        if (leftrec<=1)&&(leftsend>0),
            active = union(active,idx);
        end     
    else 
        if isempty(node(id).forward),
            if isempty(node(id).value),
                if isempty(node(id).parentmsg),
                    %node(id).forward = node(id).w(1:card(id));  %check
                    m = node(id).w*[eye(card(id));1/card(id)*ones(1,card(id))];
                    node(id).forward = m;
                    s = card(id)/node(id).w(end);
                    node(id).forwardc = (diag(m)-m'*m)/(s+1);
                else
                    %node(id).forward = forwardprop(node(id).w,node(id).parentmsg,card(node(id).parents));
                    [node(id).forward,node(id).forwardc] = forwardprop_so2(node(id).w,node(id).parentmsg,node(id).parentmsgc);
                end
            else
                %node(id).forward = zeros(1,card(id)+1);
                %node(id).forward(node(id).value+1) = 1;
                node(id).forward = zeros(1,card(id));
                node(id).forward(node(id).value+1) = 1;
                node(id).forwardc = zeros(card(id),card(id));
            end
        end
        
        node(id).sendchildren(sendid) = 1;
        idx = node(id).children(sendid);
        loc = find(node(idx).parents==id); 
        nc = length(node(id).childmsg);
        m_in = cell(nc,1);
        cov_in = cell(nc,1);
        m_in{1} = node(id).forward;
        cov_in{1} = node(id).forwardc;
        cnt = 1;
        for jj=[1:sendid-1 sendid+1:nc],
            cnt = cnt+1;
            m_in{cnt} = node(id).childmsg{jj};
            cov_in{cnt} = node(id).childmsgc{jj};
        end
        [node(idx).parentmsg{loc},node(idx).parentmsgc{loc}]= fusion_so(m_in,cov_in);    
        %node(idx).parentmsg(loc,1:card(id)+1) = approx_likefusion_all([node(id).forward; node(id).childmsg([1:sendid-1 sendid+1:end],:)]);
        node(idx).recparents(loc) = 1;
    end
    
    leftsend = sum(node(idx).sendparents==0)+sum(node(idx).sendchildren==0);
    leftrec = sum(node(idx).recchildren==0)+sum(node(idx).recparents==0);
    if (leftrec<=1)&&(leftsend>0),
        active = union(active,idx);
    end  
    
    leftsend = sum(node(id).sendparents==0)+sum(node(id).sendchildren==0);
    leftrec = sum(node(id).recparents==0)+sum(node(id).recchildren==0);
    if (leftrec==1)||(leftsend==0),
        %loc = find(active==id);
        %active(loc) = [];
        active(active==id) = [];
    end
    
end
    
    
    
for id=1:N,
    if isempty(node(id).forward),
        %node(id).forward = forwardprop(node(id).w,node(id).parentmsg,card(node(id).parents));
        [node(id).forward,node(id).forwardc] = forwardprop_so2(node(id).w,node(id).parentmsg,node(id).parentmsgc);
    end
    if isempty(node(id).back),
        %node(id).back = approx_likefusion_all(node(id).childmsg);
        [node(id).back, node(id).backc] = fusion_so(node(id).childmsg,node(id).childmsgc);
    end
    %[node(id).we,node(id).wb] = approx_likefusion_all([node(id).forward; node(id).back]);
    m_in = cell(2,1);
    cov_in = cell(2,1);
    m_in{1} = node(id).forward;
    cov_in{1} = node(id).forwardc;
    m_in{2} = node(id).back;
    cov_in{2} = node(id).backc;
    [node(id).me,node(id).mec] = fusion_so(m_in,cov_in);
end