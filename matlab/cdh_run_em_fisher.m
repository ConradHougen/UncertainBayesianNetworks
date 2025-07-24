
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

emfisherwesoinf = zeros(mxsz1, 3);
pgtemfisher = zeros(mxsz1, 1);
lemfisher = 0;

pemfisher = zeros(Nparams,NRUNS);
cemfisher = zeros(Nparams,Nparams,NRUNS);

emfisherrt = zeros(NRUNS,1);

for i = 1:NRUNS
   meta = out.meta(i);
   metafisher = fishersetup(meta,out.node(i,:));
   
   t0 = tic;
   [pemfisher(:,i),cemfisher(:,:,i)]=em_fisher(meta,out.val(:,:,i),metafisher.vala);
   emfisherrt(i) = toc(t0);
   
   meta.m = pemfisher(:,i);
   meta.cov = cemfisher(:,:,i); 
    
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
       emfisherwesoinf(1+lemfisher:lemfisher+lweso,:) = fixw(weso);

       pg = p(nzi,:,j);
       pg = pg(:);
       pgtemfisher(1+lemfisher:lemfisher+lweso,:) = pg; 
       
       lemfisher = lemfisher + lweso;
   end
        
end

avgemfisherrt = mean(emfisherrt);