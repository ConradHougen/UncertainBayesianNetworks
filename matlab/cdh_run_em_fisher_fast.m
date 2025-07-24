
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

emffwesoinf = zeros(mxsz1, 3);
pgtemff = zeros(mxsz1, 1);
lemff = 0;

pemff = zeros(Nparams,NRUNS);
cemff = zeros(Nparams,Nparams,NRUNS);

emffrt = zeros(NRUNS,1);

for i = 1:NRUNS
    fprintf(1,'EM-FF run %d\n', i);
    
    meta = out.meta(i);
    metafisher = fishersetup(meta,out.node(i,:));

    t0 = tic;
    [pemff(:,i),cemff(:,:,i)]=em_fisher_fast(meta,out.val(:,:,i),metafisher);
    emffrt(i) = toc(t0);

    meta.m = pemff(:,i);
    meta.cov = cemff(:,:,i); 

    meta = setup(out.circ, meta);

    [p,m] = inferProbfast(meta,out.val(:,:,i));

    for j = 1:Ntrain
        [pso,vso] = inferSOProbfast(meta,out.val(j,:,i));
        nzi = find(out.obsmask(j,:,i) == 0);
        pso = pso(nzi,:);
        vso = vso(nzi,:);
        pso = pso(:);
        vso = vso(:);
        weso = cdh_gen_subjective_opinions(pso,vso);

        lweso = size(weso, 1);
        emffwesoinf(1+lemff:lemff+lweso,:) = fixw(weso);

        pg = p(nzi,:,j);
        pg = pg(:);
        pgtemff(1+lemff:lemff+lweso,:) = pg; 

        lemff = lemff + lweso;
    end
        
end

avgemffrt = mean(emffrt);