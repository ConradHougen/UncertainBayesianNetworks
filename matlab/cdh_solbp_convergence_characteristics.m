% File for plotting curves of convergence of means vs. number of SOLBP
% iterations
clear;

rg = [2 3];
mxrg = rg(2);

Ntrain = 50;
Nmonte = 100;
Nnetworks = 10;
Nruns = Nmonte*Nnetworks;
SOLBPITERS = 2; % Change this to control num SOLBP iterations

node = makenetwork('stacked_diamond.file');
Nvars = length(node);

end_loc = getends(node);%the starting list of active nodes - those that have one parent or
Nends = length(end_loc);
interior_loc = setdiff((1:Nvars)',end_loc);% the list of non-active nodes
Nint = length(interior_loc);

% observed(end) and unobserved(int) variables
end_loc = [1,7];
Nends = 2;
interior_loc = [2,3,5,6,4];
Nint = 5;

% 21node net variables
% end_loc = [3, 5, 11, 10, 13, 21, 19, 16, 20, 9];
% Nends = 10;
% interior_loc = [1, 2, 4, 6, 7, 8, 12, 14, 15, 17, 18];
% Nint = 11;

% 27nodenet variables
% end_loc = [3, 5, 11, 10, 13, 21, 19, 16, 20, 9, 27, 24, 23];
% Nends = 13;
% interior_loc = [1, 2, 4, 6, 7, 8, 12, 14, 15, 17, 18, 22, 26, 25];
% Nint = 14;

obs = zeros(Nends,2);
card = zeros(1,Nvars);
dof = zeros(1,Nvars);

% Store all inferred means/variances (up to DOF) i.e. for the unobserved
% variables
pgt = zeros(Nruns*(Nvars-Nends), 1);
infpso = zeros(Nruns*(Nvars-Nends), 1);
infvso = zeros(Nruns*(Nvars-Nends), 1);
infpsolbp = zeros(Nruns*(Nvars-Nends), SOLBPITERS);
infvsolbp = zeros(Nruns*(Nvars-Nends), SOLBPITERS);
idx0 = 1;
idx1 = 1;
idx2vec = ones(1, SOLBPITERS);

for i=1:Nruns
    if mod(i-1, Nmonte)==0
        node = createBN(node,rg);
        for j = 1:Nvars
            card(j) = node(j).cardinality;
        end
        dof = card - 1;
    end
    loc = getends(node);
    [node,val] = createSLBN(node,Ntrain);

    for j=1:Nends
        value = sum(rand(1)>(1:card(end_loc(j))-1)/card(end_loc(j)));
        obs(j,:) = [end_loc(j) 1+value]; %random observation for the end nodes
        node(end_loc(j)).value = value;
    end
    
    % Create SPN for ground truth comparison
    [circ,meta] = net2circuit(node);

    % First-order inference using SPN
    p1 = inferProb(circ,meta,obs);
    
    for j = 1:length(interior_loc)
        k = interior_loc(j);
        pgt(idx0:idx0+dof(k)-1) = p1(k,1:dof(k));
        idx0 = idx0 + dof(k);
    end

    % Second-order inference using SPN
    [pso,vso] = inferSOProb(circ,meta,obs);
    for j = 1:length(interior_loc)
        k = interior_loc(j);
        infpso(idx1:idx1+dof(k)-1) = pso(k,1:dof(k));
        infvso(idx1:idx1+dof(k)-1) = vso(k,1:dof(k));
        idx1 = idx1 + dof(k); 
    end

    % Second-order inference using Loopy BP
    nodeLBP_so = inferSOloopy_v2(node);
    lennodeLBPso = length(nodeLBP_so);
    for l=1:Nvars
        card(l) = nodeLBP_so(l).cardinality;
        nodeLBP_so(l).val = [];
        nodeLBP_so(l).childmsg = cell(length(nodeLBP_so(l).children),1);
        nodeLBP_so(l).childmsgc = cell(length(nodeLBP_so(l).children),1); 
        nodeLBP_so(l).parentmsg = cell(length(nodeLBP_so(l).parents),1);
        nodeLBP_so(l).parentmsgc = cell(length(nodeLBP_so(l).parents),1);
        nodeLBP_so(l).forward = [];
        nodeLBP_so(l).forwardc = [];
        nodeLBP_so(l).back = [];
        nodeLBP_so(l).backc = [];
        nodeLBP_so(l).delta_me = 0;
        nodeLBP_so(l).me = zeros(1,card(l));
        
        % Init the childmsg and parentmsg fields
        % Init the childmsgc and parentmsgc fields
        nc = length(nodeLBP_so(l).children);
        for j = 1:nc
            nodeLBP_so(l).childmsg{j} = ones(1,card(l))/card(l);
            nodeLBP_so(l).childmsgc{j} = zeros(card(l),card(l));
        end
        np = length(nodeLBP_so(l).parents);
        for j = 1:np
            pcrd = card(nodeLBP_so(l).parents(j));
            nodeLBP_so(l).parentmsg{j} = ones(1,pcrd);
            nodeLBP_so(l).parentmsgc{j} = zeros(pcrd,pcrd);
        end
        
    end
    
    for l=1:size(obs,1)
        nodeLBP_so(obs(l,1)).val = obs(l,2);
    end

    % BEGIN SOLBP ITERATIONS
    niter = 0;
    while niter < SOLBPITERS
        active = [1:Nvars];
        id = active(randi(length(active)));
        
        % Loop through nodes in BN, sending the parent/child messages
        while ~isempty(active)
            % Send parent messages
            if isempty(nodeLBP_so(id).val)
                if isempty(nodeLBP_so(id).childmsg)
                    crd = card(id);
                    nodeLBP_so(id).back = ones(1,crd)./crd;
                    nodeLBP_so(id).backc = zeros(crd,crd);
                else
                    [nodeLBP_so(id).back, nodeLBP_so(id).backc] = fusion_so(nodeLBP_so(id).childmsg,nodeLBP_so(id).childmsgc);
                    % DEBUG
    %                 if id == 3
    %                     fprintf(1,"node(3).back=%1.4f,$1.4f\n",node(id).back(1),node(id).back(2));
    %                 end
                end
            else
                crd = card(id);
                nodeLBP_so(id).back = zeros(1,crd);
                nodeLBP_so(id).back(nodeLBP_so(id).val) = 1;
                %node(id).back(node(id).val+1) = 1;
                nodeLBP_so(id).backc = zeros(crd,crd);
            end
            for sendid = 1:length(nodeLBP_so(id).parents)
                idx = nodeLBP_so(id).parents(sendid);
                loc = find(nodeLBP_so(idx).children == id);
                [nodeLBP_so(idx).childmsg{loc},nodeLBP_so(idx).childmsgc{loc}] = backprop_so(nodeLBP_so(id).w,nodeLBP_so(id).parentmsg,nodeLBP_so(id).parentmsgc,card(nodeLBP_so(id).parents),sendid,nodeLBP_so(id).back,nodeLBP_so(id).backc);
            end
    
            % Send child messages
            if isempty(nodeLBP_so(id).val)
                [nodeLBP_so(id).forward,nodeLBP_so(id).forwardc] = forwardprop_so(nodeLBP_so(id).w,nodeLBP_so(id).parentmsg,nodeLBP_so(id).parentmsgc);
                %DEBUG
    %             if id == 3
    %                 fprintf(1,"node(3).forward=%1.4f,%1.4f\n",node(id).forward(1),node(id).forward(2));
    %             end
            else
    %             m = node(id).w*[eye(card(id));1/card(id)*ones(1,card(id))];
    %             node(id).forward = m;
    %             s = card(id)/node(id).w(end);
    %             node(id).forwardc = (diag(m)-m'*m)/(s+1);
                nodeLBP_so(id).forward = zeros(1,card(id));
                nodeLBP_so(id).forward(nodeLBP_so(id).val) = 1;
                nodeLBP_so(id).forwardc = zeros(card(id),card(id));
            end
            for sendid = 1:length(nodeLBP_so(id).children)
                idx = nodeLBP_so(id).children(sendid);
                loc = find(nodeLBP_so(idx).parents==id);
                nc = length(nodeLBP_so(id).childmsg);
                m_in = cell(nc,1);
                cov_in = cell(nc,1);
                m_in{1} = nodeLBP_so(id).forward;
                cov_in{1} = nodeLBP_so(id).forwardc;
                cnt = 1;
                for jj=[1:sendid-1 sendid+1:nc]
                    cnt = cnt+1;
                    m_in{cnt} = nodeLBP_so(id).childmsg{jj};
                    cov_in{cnt} = nodeLBP_so(id).childmsgc{jj};
                end
                [nodeLBP_so(idx).parentmsg{loc},nodeLBP_so(idx).parentmsgc{loc}]= fusion_so(m_in,cov_in);
            end
    
            % Remove id from active list and get the next id
            loc = active==id;
            active(loc) = [];
            if ~isempty(active)
                id = active(randi(length(active)));
            end
        end
        maxdelta = 0;
        for id=1:Nvars
            if isempty(nodeLBP_so(id).forward)
                [nodeLBP_so(id).forward,nodeLBP_so(id).forwardc] = forwardprop_so(nodeLBP_so(id).w,nodeLBP_so(id).parentmsg,nodeLBP_so(id).parentmsgc);
            end
            if isempty(nodeLBP_so(id).back)
                [nodeLBP_so(id).back, nodeLBP_so(id).backc] = fusion_so(nodeLBP_so(id).childmsg,nodeLBP_so(id).childmsgc);
            end
            prev_me = nodeLBP_so(id).me;
            
            m_in = cell(2,1);
            cov_in = cell(2,1);
            m_in{1} = nodeLBP_so(id).forward;
            cov_in{1} = nodeLBP_so(id).forwardc;
            m_in{2} = nodeLBP_so(id).back;
            cov_in{2} = nodeLBP_so(id).backc;
            [nodeLBP_so(id).me,nodeLBP_so(id).mec] = fusion_so(m_in,cov_in);
            
            nodeLBP_so(id).delta_me = max(abs(prev_me - nodeLBP_so(id).me));
            
            maxdelta = max(maxdelta, nodeLBP_so(id).delta_me);
        end
        
        idx2 = idx2vec(niter+1);
        for j = 1:length(interior_loc)
            k = interior_loc(j);
            infpsolbp(idx2:idx2+dof(k)-1, niter+1) = nodeLBP_so(k).me(1:dof(k));
            mecdiag = diag(nodeLBP_so(k).mec);
            infvsolbp(idx2:idx2+dof(k)-1, niter+1) = mecdiag(1:dof(k));
            idx2 = idx2 + dof(k); 
        end
        idx2vec(niter+1) = idx2;

        niter = niter + 1;
    end
    % End of SOLBP loop
    
end

weso = cdh_gen_subjective_opinions(infpso,infvso);
wesolbp = cdh_gen_subjective_opinions(infpsolbp(:,end),infvsolbp(:,end));

% Compare inferred values from SPN and LBP methods
figure;
scatter(infpso, infpsolbp(:,end), 'MarkerEdgeColor','b');
title('Inferred Means: ');
xlabel('Mean Value (SPN Method)');
ylabel('Mean Value (SOLBP)');

figure;
scatter(infvso, infvsolbp(:,end), 'MarkerEdgeColor','b');
title('Inferred Variances: ');
xlabel('Variance Value (SPN Method)');
ylabel('Variance Value (SOLBP)');

% Figures showing avg errors (means and vars) vs. number of SOLBP
% iterations
avgabserrspso = zeros(1, SOLBPITERS);
avgabserrsvso = zeros(1, SOLBPITERS);
for r = 1:SOLBPITERS
    avgabserrspso(r) = mean(abs(infpsolbp(:,r)-infpso));
    avgabserrsvso(r) = mean(abs(infvsolbp(:,r)-infvso));
end
figure;
scatter(1:SOLBPITERS, avgabserrspso, 'MarkerEdgeColor','b');
title("Avg Abs Errs of Means vs. SOLBP Iteration");
xlabel("SOLBP Iteration");
ylabel("Avg Abs Errs of Means");

figure;
scatter(1:SOLBPITERS, avgabserrsvso, 'MarkerEdgeColor','b');
title("Avg Abs Errs of Vars vs. SOLBP Iteration");
xlabel("SOLBP Iteration");
ylabel("Avg Abs Errs of Vars");

% Create scatter plot comparing absolute errors of means to absolute error
% of variances (test for correlation)
absdiffspso = abs(infpsolbp(:,end) - infpso);
absdiffsvso = abs(infvsolbp(:,end) - infvso);
randabsdiffspso = absdiffspso(randperm(length(absdiffspso)));
randabsdiffsvso = absdiffsvso(randperm(length(absdiffsvso)));
corr1 = corr(absdiffspso, absdiffsvso);
corr2 = corr(randabsdiffspso, randabsdiffsvso);
figure;
scatter(absdiffspso, absdiffsvso, 'MarkerEdgeColor','b');
title("Abs Var Errors vs. Abs Mean Errors (Corr "+corr1+" > "+corr2+")");
xlabel('Absolute Mean Error');
ylabel('Absolute Variance Error');

[pa,pp] = perfbound_beta_new2(pgt, weso);
[pa2,pp2] = perfbound_beta_new2(pgt, wesolbp);
figure; 
plot(pp, pa, 'b', pp2, pa2, 'r+', pp2, pp2, 'k--');
title('DeCBoD: ');
xlabel('Desired Confidence');
ylabel('Actual Confidence');
legend('SPN Method', 'SOLBP', 'Reference Line', 'Location', 'northwest');