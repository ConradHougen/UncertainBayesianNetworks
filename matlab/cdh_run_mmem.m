
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

mmemwesoinf = zeros(mxsz1, 3);
pgtmmem = zeros(mxsz1, 1);
lmmem = 0;

pmmem = zeros(Nparams,NRUNS);
cmmem = zeros(Nparams,Nparams,NRUNS);

for i = 1:NRUNS
   meta = out.meta(i);
   [pmmem(:,i),cmmem(:,:,i),~]=mmem(out.circ,meta,out.val(:,:,i)); 
   
   meta.m = pmmem(:,i);
   meta.cov = cmmem(:,:,i); 

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
       mmemwesoinf(1+lmmem:lmmem+lweso,:) = fixw(weso);
       
       pg = p(nzi,:,j);
       pg = pg(:);
       pgtmmem(1+lmmem:lmmem+lweso,:) = pg;
       
       lmmem = lmmem + lweso;
   end
end

