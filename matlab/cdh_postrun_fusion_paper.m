% Postrun script for case of 2-variable binary Bayes Net. This script
% computes the Fusion paper results: learnposterior computes and stores the
% means and covariances via counting, storing them in minct and rinct.
% Using these means and covariances, we would like to perform second order
% inference to generate the subjective opinions. We then plot the
% desired confidence bound divergence curves, where perfbound is computed
% between the pgtinfer from inferProb and the subjective opinions 

% Mask for Independent Fisher (orange curve in Fusion)
covmsk = zeros(Nparams,Nparams);
for i = 1:size(covmsk,1)/2
    covmsk(1+(i-1)*2:i*2,1+(i-1)*2:i*2) = ones(2,2);
end

% minct and rinct are in terms of free params only -> need to expand these
% back out to all params in order to use inferSOProb

fusionweso = zeros(linfactive, 3);
zfusionweso = zeros(linfactive, 3);
lfusionweso = 0;
mincall = zeros(Nparams*NRUNS,1);
drincall = zeros(Nparams*NRUNS,1);
% Populate fusionweso using the full covariance rinct
for i = 1:NRUNS
    
    if mod(i, 10) == 0
        fprintf(1,'Starting Run %d\n',i);
    end
    
    % First, expand minct and rinct into Nparams from Nfreeparams
    minc = zeros(Nparams, 1);
    rinc = zeros(Nparams, Nparams);
    for ii = 1:Nfreeparams
        minc(2*(ii-1)+1) = minct(ii,i);
        minc(2*(ii-1)+2) = 1 - minct(ii,i);
      
        for jj = 1:Nfreeparams
            rinc(2*(ii-1)+1,2*(jj-1)+1) = rinct(ii,jj,i);
            rinc(2*(ii-1)+1,2*(jj-1)+2) = -rinct(ii,jj,i);
            rinc(2*(ii-1)+2,2*(jj-1)+1) = -rinct(ii,jj,i);
            rinc(2*(ii-1)+2,2*(jj-1)+2) = rinct(ii,jj,i);
        end
    end
    
    drinc = diag(rinc); %Only the variances
    mincall((i-1)*Nparams+1:i*Nparams) = minc;
    drincall((i-1)*Nparams+1:i*Nparams) = drinc;
    
    meta2 = metaarr{floor(i/NDATASETS) + 1};
    meta2.m = minc;
    for j = 1:Ntrain
        v = val(j,:,i);
        obs = [(1:Nvars)' v'];
        obsidx = find(obs(:,2) ~= 0);
        obs = obs(obsidx,:);
        nobsidx = setdiff(1:Nvars, obsidx);
        
  
        % Full covariance
        meta2.cov = rinc;
        [pso, vso] = inferSOProb(circ,meta2,obs);
        pso = pso(nobsidx, :);
        vso = vso(nobsidx, :);
        pso = pso(:);
        vso = vso(:);
        weso = cdh_gen_subjective_opinions(pso, vso);
        lweso = size(weso, 1);
        fusionweso(lfusionweso+1:lfusionweso+lweso,:) = weso;
        
        % Zero off-diagonal covariance
        meta2.cov = rinc .* covmsk;
        [pso, vso] = inferSOProb(circ,meta2,obs);
        pso = pso(nobsidx, :);
        vso = vso(nobsidx, :);
        pso = pso(:);
        vso = vso(:);
        weso = cdh_gen_subjective_opinions(pso, vso);
        zfusionweso(lfusionweso+1:lfusionweso+lweso,:) = weso;
        
        lfusionweso = lfusionweso + lweso;
    end
    
end

xline = 0:0.01:1;

% DeCBod prior to inferSOProb
weso = cdh_gen_subjective_opinions(mincall, drincall);
[paFnoInfSO,ppFnoInfSO] = perfbound_beta_new2(pgt,weso);
figure;
plot(paFnoInfSO, ppFnoInfSO);
hold on;
plot(xline, xline);
title('DeCBoD (Fusion mean/var only)');
xlabel('pa');
ylabel('pp');
        

% Threshold entries near zero
fusionweso(abs(fusionweso) < ZERO) = 0;
zfusionweso(abs(zfusionweso) < ZERO) = 0;

% Desired confidence bound divergence for full covariance
[paF,ppF] = perfbound_beta_new2(pgtinfer(1:lfusionweso), fusionweso(1:lfusionweso,:));
[paZF,ppZF] = perfbound_beta_new2(pgtinfer(1:lfusionweso), zfusionweso(1:lfusionweso,:));
figure;
plot(paF,ppF);
hold on;
plot(xline,xline,'--');
plot(paZF,ppZF);
title('DeCBoD Fusion Paper (3-node chain)');
xlabel('pa');
ylabel('pp');

