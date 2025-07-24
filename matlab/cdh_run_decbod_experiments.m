
clear;
% 1. Set up experiment parameters
rg = [2 2];
netfile = '3node_1_1_1.file';
Ntrain = 120; % Number of training examples in val matrix
PARTIALOBS = true; % Only make partial observations?
FRACOBS = 0.5; % What fraction [0,1] of observations should we keep?
NGTBNS = 1;
NDATASETS = 1000;
NRUNS = NGTBNS*NDATASETS;

% 2. Generate ground truth bayes nets and observations
out = cdh_generate_gtbns_vals(rg,netfile,NRUNS,Ntrain,FRACOBS);
Nparams = length(out.meta(1).varmap);

% 3. Run Methods
xline = 0:0.01:1;

% Method 1: Method of Moments
cdh_run_mm;
[pamm,ppmm] = perfbound_beta_new2(pgtm(1:lmm), mmwesoinf(1:lmm,:));
rmsemm = cdh_rms_error(pgtm,mmwesoinf);
avgdecbod_mm = sum(abs(pamm - ppmm))/length(pamm);

% Method 2: Expectation Maximization
% cdh_run_em;
% [paem,ppem] = perfbound_beta_new2(pgte(1:lem), emwesoinf(1:lem,:));
% figure;
% plot(ppem,paem);
% hold on;
% xline = 0:0.01:1;
% plot(xline,xline,'--');
% title('DeCBoD (EM Fast)');
% xlabel('pp');
% ylabel('pa');
% rmseem = cdh_rms_error(pgte,emwesoinf);

% Method 3: EM with Gaussian Approximation
cdh_run_em_fast_map;
[paemmap,ppemmap] = perfbound_beta_new2(pgtemmap(1:lemmap), emmapwesoinf(1:lemmap,:));
rmseemmap = cdh_rms_error(pgtemmap,emmapwesoinf);
avgdecbod_emmap = sum(abs(paemmap - ppemmap))/length(paemmap);

% Method 4: EM Fisher
% cdh_run_em_fisher;
% [paemfisher,ppemfisher] = perfbound_beta_new2(pgtemfisher(1:lemfisher), emfisherwesoinf(1:lemfisher,:));
% figure;
% plot(ppemfisher,paemfisher);
% hold on;
% xline = 0:0.01:1;
% plot(xline,xline,'--');
% title('DeCBoD (EM-Fisher)');
% xlabel('pp');
% ylabel('pa');
% rmseemfisher = cdh_rms_error(pgtemfisher,emfisherwesoinf);

% Method 5: EM Fisher Fast
cdh_run_em_fisher_fast;
[paemff,ppemff] = perfbound_beta_new2(pgtemff(1:lemff), emffwesoinf(1:lemff,:));
rmseemff = cdh_rms_error(pgtemff,emffwesoinf);
avgdecbod_emff = sum(abs(paemff - ppemff))/length(paemff);

% langeparams.MAXITER = 15000;
% langeparams.LANGEMAXSAMP = 6000;
% langeparams.LANGEBURNIN = 2000;
% langeparams.LANGEALPHATHRESH = 0.8;
% langeparams.Ntrain = Ntrain;
% langeparams.Nminibatch = 1;
% langeparams.fpidx = out.fpidx;
% langeparams.ifpidx = out.ifpidx;
% langeparams.circ = out.circ;
% cdh_run_ldmh;
% [paldmh,ppldmh] = perfbound_beta_new2(pgtldmh(1:lldmh), ldmhwesoinf(1:lldmh,:));
% figure;
% plot(ppldmh,paldmh);
% hold on;
% plot(xline,xline,'--');
% title('DeCBoD (LDMH)');
% xlabel('pp');
% ylabel('pa');


% figure;
% newcolors = [0.83 0.14 0.14
%              1.00 0.54 0.00
%              0.47 0.25 0.80
%              0.25 0.80 0.54];
% colororder(newcolors)
% plot(xline,xline,'--k');
% hold on;
% title('DeCBoD: 9-node Inverted Trees');
% xlabel('Desired Confidence');
% ylabel('Actual Confidence');
% 
% plot(ppmm,pamm,'-','LineWidth',1.2);
% plot(ppemmap,paemmap,'-','LineWidth',1.2);
% plot(ppemff,paemff,'-','LineWidth',1.2);
% 
% legend('y=x','BMM','EM-GA','EM-Fisher','Location','northwest');