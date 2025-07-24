
mxsz1 = ceil(Ntrain*rg(2)*out.Nvars*NRUNS*(1-FRACOBS));

emwesoinf = zeros(mxsz1, 3);
pgte = zeros(mxsz1, 1);
lem = 0;

pem = zeros(Nparams,NRUNS);
cem = zeros(Nparams,Nparams,NRUNS);

emrt = zeros(NRUNS,1);

for i = 1:NRUNS
   meta = out.meta(i);
   
   t0 = tic;
   [pem(:,i),cem(:,:,i)]=em_fastv2(meta,out.val(:,:,i));
   emrt(i) = toc(t0);
   
   meta.m = pem(:,i);
   meta.cov = cem(:,:,i); 
    
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
       emwesoinf(1+lem:lem+lweso,:) = fixw(weso);

       pg = p(nzi,:,j);
       pg = pg(:);
       pgte(1+lem:lem+lweso,:) = pg; 
       
       lem = lem + lweso;
   end
        
end

avgemrt = mean(emrt);