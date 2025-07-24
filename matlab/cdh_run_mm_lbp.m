% Test file for testing Loopy Belief Prop for second-order learning using Method of Moments
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

mmwesoinf = zeros(mxsz1, 3);
pgtm = zeros(mxsz1, 1);
lmm = 0;
pmm = zeros(Nparams,NRUNS);
cmm = zeros(Nparams,Nparams,NRUNS);
mmrt = zeros(NRUNS,1);

for i = 1:NRUNS
    fprintf(1,'MM run %d\n', i);
    meta = out.meta(i);
    node = out.node(i,:);

    t0 = tic;
    [pmm(:,i),cmm(:,:,i),~]= mm_lbp(meta,node,out.val(:,:,i));
    mmrt(i) = toc(t0);

    meta.m = pmm(:,i);
    meta.cov = cmm(:,:,i);

    [p,m] = inferProbfast(meta,out.val(:,:,i));

    for j = 1:Ntrain
        [pso,vso] = inferSOProbfast(meta,out.val(j,:,i));
        nzi = find(out.obsmask(j,:,i) == 0);
        pso = pso(nzi,:);
        vso = vso(nzi,:);
        pso = pso(:);
        vso = vso(:);
        weso = cdh_gen_subjective_opinions(pso,vso);

        lweso = size(weso,1);       
        mmwesoinf(1+lmm:lmm+lweso,:) = fixw(weso);

        pg = p(nzi,:,j);
        pg = pg(:);
        pgtm(1+lmm:lmm+lweso,:) = pg;

        lmm = lmm + lweso;
    end

end

avgmmrt = mean(mmrt);
