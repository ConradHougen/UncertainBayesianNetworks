clear;
for jj = 1:20

    % File for running tests related to SO inferencing with LBP
    close all;

    rg = [2 2];
    mxrg = rg(2);

    Ntrain = 10000;
    Nmonte = 50;
    Nnetworks = 10;
    Nruns = Nmonte*Nnetworks;

    node = makenetwork('21nodenet.file');
    Nvars = length(node);

    end_loc = getends(node);%the starting list of active nodes - those that have one parent or
    Nends = length(end_loc);
    interior_loc = setdiff((1:Nvars)',end_loc);% the list of non-active nodes
    Nint = length(interior_loc);

    % observed(end) and unobserved(int) variables
    % end_loc = [2,3];
    % Nends = 2;
    % interior_loc = [1,4,5,6];
    % Nint = 4;

    % 21node net variables
    end_loc = [3, 5, 11, 10, 13, 21, 19, 16, 20, 9];
    Nends = 10;
    interior_loc = [1, 2, 4, 6, 7, 8, 12, 14, 15, 17, 18];
    Nint = 11;

    % 27nodenet variables
    % end_loc = [3, 5, 11, 10, 13, 21, 19, 16, 20, 9, 27, 24, 23];
    % Nends = 13;
    % interior_loc = [1, 2, 4, 6, 7, 8, 12, 14, 15, 17, 18, 22, 26, 25];
    % Nint = 14;

    obs = zeros(Nends,2);

    %pgt = zeros(Nruns*Nint*mxrg,1);
    %pgtLBP = zeros(Nruns*Nint*mxrg,1);
    weso = zeros(Nruns*Nint*mxrg,3);
    wesolbp = zeros(Nruns*Nint*mxrg,3);

    tm0 = zeros(Nruns,1);
    tm_slbn = zeros(Nruns,1);

    tm_so = zeros(Nruns,1);
    tm_fo = zeros(Nruns,1);
    tm_folbp = zeros(Nruns,1);
    tm_solbp = zeros(Nruns,1);

    card = zeros(1,Nvars);
    dof = zeros(1,Nvars);

    % Store all inferred means/variances (up to DOF) i.e. for the unobserved
    % variables
    pgt = zeros(Nruns*(Nvars-Nends), 1);
    infpso = zeros(Nruns*(Nvars-Nends), 1);
    infvso = zeros(Nruns*(Nvars-Nends), 1);
    infpsolbp = zeros(Nruns*(Nvars-Nends), 1);
    infvsolbp = zeros(Nruns*(Nvars-Nends), 1);
    idx0 = 1;
    idx1 = 1;
    idx2 = 1;

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

        for j=1:Nends,
            value = sum(rand(1)>(1:card(end_loc(j))-1)/card(end_loc(j)));
            obs(j,:) = [end_loc(j) 1+value]; %random observation for the end nodes
            node(end_loc(j)).value = value;
        end

        % Create SPN for ground truth comparison
        tic;
        [circ,meta] = net2circuit(node);
        tm0(i) = toc;

        % First-order inference using SPN
        tic;
        p1 = inferProb(circ,meta,obs);
        tm_fo(i) = toc;

        for j = 1:length(interior_loc)
            k = interior_loc(j);
            pgt(idx0:idx0+dof(k)-1) = p1(k,1:dof(k));
            idx0 = idx0 + dof(k);
        end

        % Second-order inference using SPN
        tic;
        [pso,vso] = inferSOProb(circ,meta,obs);
        tm_so(i) = toc;
        for j = 1:length(interior_loc)
            k = interior_loc(j);
            infpso(idx1:idx1+dof(k)-1) = pso(k,1:dof(k));
            infvso(idx1:idx1+dof(k)-1) = vso(k,1:dof(k));
            idx1 = idx1 + dof(k); 
        end


    %     pso = pso';
    %     vso = vso';
    %     pso = pso(:);
    %     vso = vso(:);
    %     weso = cdh_gen_subjective_opinions(pso,vso);

        % First-order inference using Loopy BP
        tic;
        nodeLBP_fo = prob_infer_loopy(node,obs);
        tm_folbp(i) = toc;

        % Second-order inference using Loopy BP
        tic;
        nodeLBP_so = inferSOloopy(node,obs,jj);
        tm_solbp(i) = toc;
        fprintf(1, "SOLBP with niter=%d, trial=%d\n", jj, i);
        sz = 0;
        for k = 1:Nvars
            sz = sz + card(k);
        end
        psoLBP = zeros(sz, 1);
        vsoLBP = zeros(sz, 1);
        next = 1;
        for k = 1:Nvars
            psoLBP(next:next+card(k)-1) = nodeLBP_so(k).me';
            vsoLBP(next:next+card(k)-1) = diag(nodeLBP_so(k).mec);
            next = next + card(k);
        end

        for j = 1:length(interior_loc)
            k = interior_loc(j);
            infpsolbp(idx2:idx2+dof(k)-1) = nodeLBP_so(k).me(1:dof(k));
            mecdiag = diag(nodeLBP_so(k).mec);
            infvsolbp(idx2:idx2+dof(k)-1) = mecdiag(1:dof(k));
            idx2 = idx2 + dof(k); 
        end

    %     wesolbp = cdh_gen_subjective_opinions(psoLBP,vsoLBP);

    end

    weso = cdh_gen_subjective_opinions(infpso,infvso);
    wesolbp = cdh_gen_subjective_opinions(infpsolbp,infvsolbp);

    % Compare inferred values from SPN and LBP methods
    gcf = figure;
    scatter(infpso, infpsolbp, 'MarkerEdgeColor','b');
    title('Inferred Means: 21-Node BN');
    xlabel('Mean Value (SPN Method)');
    ylabel('Mean Value (SOLBP)');
    fstr = ".\twentyonenodenetv2\infmeans_niter_" + string(jj) + ".tif";
    saveas(gcf, fstr);
    fstr = ".\twentyonenodenetv2\infmeans_niter_" + string(jj) + ".fig";
    saveas(gcf, fstr);

    gcf = figure;
    scatter(infvso, infvsolbp, 'MarkerEdgeColor','b');
    title('Inferred Variances: 21-Node BN');
    xlabel('Variance Value (SPN Method)');
    ylabel('Variance Value (SOLBP)');
    fstr = ".\twentyonenodenetv2\infvars_niter_" + string(jj) + ".tif";
    saveas(gcf, fstr);
    fstr = ".\twentyonenodenetv2\infvars_niter_" + string(jj) + ".fig";
    saveas(gcf, fstr);
    

    % Create scatter plot comparing absolute errors of means to absolute error
    % of variances (test for correlation)
    % absdiffspso = abs(infpsolbp - infpso);
    % absdiffsvso = abs(infvsolbp - infvso);
    % randabsdiffspso = absdiffspso(randperm(length(absdiffspso)));
    % randabsdiffsvso = absdiffsvso(randperm(length(absdiffsvso)));
    % corr1 = corr(absdiffspso, absdiffsvso);
    % corr2 = corr(randabsdiffspso, randabsdiffsvso);
    % figure;
    % scatter(absdiffspso, absdiffsvso, 'MarkerEdgeColor','b');
    % title("Abs Var Errors vs. Abs Mean Errors (Corr "+corr1+" > "+corr2+")");
    % xlabel('Absolute Mean Error');
    % ylabel('Absolute Variance Error');

    [pa,pp] = perfbound_beta_new2(pgt, weso);
    [pa2,pp2] = perfbound_beta_new2(pgt, wesolbp);
    gcf = figure; 
    plot(pp, pa, 'b', pp2, pa2, 'r+', pp2, pp2, 'k--');
    title('DeCBoD: 21-Node BN');
    xlabel('Desired Confidence');
    ylabel('Actual Confidence');
    legend('SPN Method', 'SOLBP', 'Reference Line', 'Location', 'northwest');
    fstr = ".\twentyonenodenetv2\decbod_niter_" + string(jj) + ".tif";
    saveas(gcf, fstr);
    fstr = ".\twentyonenodenetv2\decbod_niter_" + string(jj) + ".fig";
    saveas(gcf, fstr);
    
end

