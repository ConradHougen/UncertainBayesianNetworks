
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

emmapwesoinf = zeros(mxsz1, 3);
pgtemmap = zeros(mxsz1, 1);
lemmap = 0;

pemmap = zeros(Nparams,NRUNS);
cemmap = zeros(Nparams,Nparams,NRUNS);

emmaprt = zeros(NRUNS,1);

for i = 1:NRUNS
    fprintf(1,'EM-GA run %d\n', i);
    
    meta = out.meta(i);

    t0 = tic;
    [pemmap(:,i),cemmap(:,:,i)]=em_fast_map(meta,out.val(:,:,i));
    emmaprt(i) = toc(t0);

    meta.m = pemmap(:,i);
    meta.cov = cemmap(:,:,i); 

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
        emmapwesoinf(1+lemmap:lemmap+lweso,:) = fixw(weso);

        pg = p(nzi,:,j);
        pg = pg(:);
        pgtemmap(1+lemmap:lemmap+lweso,:) = pg; 

        lemmap = lemmap + lweso;
    end
        
end

avgemmaprt = mean(emmaprt);