
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

ldmhwesoinf = zeros(mxsz1, 3);
pgtldmh = zeros(mxsz1, 1);
lldmh = 0;

pldmh = zeros(Nparams,NRUNS);
cldmh = zeros(Nparams,Nparams,NRUNS);

for i = 1:NRUNS
    meta = out.meta(i);
    
    fprintf(1, '\n---STARTING RUN %d---\n', i);
    [pldmh(:,i),cldmh(:,:,i)] = ldmh(out,meta,out.val(:,:,i),langeparams,i);
    
    meta.m = pldmh(:,i);
    meta.cov = cldmh(:,:,i);
    
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
       
       lweso = size(weso,1);
       ldmhwesoinf(1+lldmh:lldmh+lweso,:) = fixw(weso);
       
       pg = p(nzi,:,j);
       pg = pg(:);
       pgtldmh(1+lldmh:lldmh+lweso,:) = pg;
       
       lldmh = lldmh + lweso;
    end
end