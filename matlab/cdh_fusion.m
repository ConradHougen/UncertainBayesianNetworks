% Main file

% Flags/Variable Definitions
rg = [2 2]; % Cardinality range for variables
Nminibatch = 50; % Number of samples per gradient step (stochastic vs. batch GD) 
Ntrain = 500;  % Number of training examples (must be multiple of Nminibatch)
PARTIALOBS = true; % Only make partial observations?
FRACOBS = 0.75; % What fraction [0,1] of observations should we keep?
MAXITER = 15000; % Absolute maximum number of updates to theta
STOP = 3; % How many steady iterates do we need to stop iterating?
CONVERGE = 1e-5; % Threshold for "convergence" in iterates of theta vector
ZERO = 1e-10; % Approximation to 0 for numerical stability purposes
LANGEFLAG = true; % Run Langevin dynamics?
LANGEMAXSAMP = 6000; % Max number of samples under Langevin dynamics phase
LANGEBURNIN = 2000; % Number of Langevin samples to discard when inferencing
LANGEALPHATHRESH = 0.8; % Param for determining entrance to Langevin phase
NGTBNS = 1;
NDATASETS = 200;
NRUNS = NGTBNS*NDATASETS; % Number of times to run Monte Carlo simulation (>= 1)

% Create Bayes Net from net file
node = makenetwork('net2node.file');

% Number of variables in the Bayes Net
Nvars = length(node);

% Create new ground truth Bayes Net
node = createBN(node, rg);
% Translate BN into sum-product network (circ/meta)
[circ, meta] = cdh_net2circuit(node);
meta2 = meta; % meta2 can be modified
% Number of params is equal to length of meta.varmap
Nparams = length(meta.varmap);
Nnodes = length(circ); % Number of nodes in sum-product network
% Determine locations of free parameters
[fpidx, ifpidx] = cdh_find_free_params(node, meta);
Nfreeparams = length(fpidx);

% inferProb gives the inferred probabilities (expected value)
pgtinfer = zeros(ceil(Ntrain*rg(2)*Nvars*NRUNS*(1-FRACOBS)*1.2), 1);
linfactive = 0;

% Save observations from each run
val = zeros(Ntrain, Nvars, NRUNS);
obsmask = true(size(val));

% For each of the Monte Carlo runs, compute means and variances. On each MC
% run, we generate a new set of observations
pgt = [];

% Store mcomp, rcomp, minc, rinc from each run
mcompt = zeros(Nfreeparams,1,NRUNS);
rcompt = zeros(Nfreeparams,Nfreeparams,NRUNS);
minct = zeros(Nfreeparams,1,NRUNS);
rinct = zeros(Nfreeparams,Nfreeparams,NRUNS);

% Save the meta from each run
metaarr = {};
metaarr{1} = meta;
for mcrun = 1:NRUNS

    fprintf(1, '\n---STARTING RUN %d---\n', mcrun);
    
    % Every K runs, create a new ground truth BN
    if mod(mcrun, NDATASETS) == 0
        fprintf(1, '++++RUN %d: GENERATING NEW GT BN\n', mcrun);
        % Create new ground truth Bayes Net
        node = createBN(node, rg);
        % Translate BN into sum-product network (circ/meta)
        [circ, metaarr{mcrun/NDATASETS + 1}] = cdh_net2circuit(node);
        % Number of params is equal to length of meta.varmap
        Nparams = length(metaarr{mcrun/NDATASETS + 1}.varmap);
        Nnodes = length(circ); % Number of nodes in sum-product network
        % Determine locations of free parameters
        [fpidx, ifpidx] = cdh_find_free_params(node, metaarr{mcrun/NDATASETS + 1});
        Nfreeparams = length(fpidx);
    end
    
    % Get ground truth probabilities from GT BN
    for i = 1:Nvars
        pgt = [pgt; node(i).p(:)];
    end

    % Generate complete set of Ntrain observations
    val(:,:,mcrun) = cdh_sim_observations(node, Ntrain);

    % Some of the nodes are unobserved (create mask for which vars are observed
    % over the Ntrain samples).
    if PARTIALOBS == true
        obsmask(:,:,mcrun) = rand(size(val(:,:,mcrun))) < FRACOBS;
    end

    % Now that we have all observations (including unobserved vars), infer mean
    % and covariance for both complete and incomplete observations
    % RETURNS 0 in cases other than 2-variable Bayes Net case! i.e. doesn't
    % work for those cases, but will not crash
    [mcompt(:,:,mcrun), rcompt(:,:,mcrun)] = learnposterior(val(:,:,mcrun));
    [minct(:,:,mcrun), rinct(:,:,mcrun)] = learnposterior(val(:,:,mcrun) .* obsmask(:,:,mcrun));

    % Apply obsmask to "remove" observations
    val(:,:,mcrun) = val(:,:,mcrun) .* obsmask(:,:,mcrun);
    
    for i = 1:Ntrain
        v = val(i,:,mcrun);
        obs = [(1:Nvars)' v'];
        obsidx = find(obs(:,2) ~= 0);
        obs = obs(obsidx,:);
        nobsidx = setdiff(1:Nvars, obsidx);
    
        % IMPORTANT: Use the original meta from net2circuit here
        ipgt = inferProb(circ,metaarr{floor(mcrun/NDATASETS) + 1},obs);
        ipgt = ipgt(nobsidx, :);
        ipgt = ipgt(:);
        lplus = length(ipgt);
        pgtinfer(linfactive+1:linfactive+lplus) = ipgt;
        
        linfactive = linfactive + lplus; 
    end
   
end



    