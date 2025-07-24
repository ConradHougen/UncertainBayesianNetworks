% Test performance bound using bernoulli variable Bayes networks

% Flags/Variable Definitions
rg = [2 2]; % Cardinality range for variables
Nminibatch = 50; % Number of samples per gradient step (stochastic vs. batch GD) 
Ntrain = 50;  % Number of training examples (must be multiple of Nminibatch)
PARTIALOBS = true; % Only make partial observations?
FRACOBS = 0.8; % What fraction [0,1] of observations should we keep?
MAXITER = 10000; % Absolute maximum number of updates to theta
STOP = 3; % How many steady iterates do we need to stop iterating?
CONVERGE = 1e-5; % Threshold for "convergence" in iterates of theta vector
ZERO = 1e-10; % Approximation to 0 for numerical stability purposes
LANGEFLAG = true; % Run Langevin dynamics?
LANGEMAXSAMP = 4000; % Max number of samples under Langevin dynamics phase
LANGEBURNIN = 1000; % Number of Langevin samples to discard when inferencing
LANGEALPHATHRESH = 0.8; % Param for determining entrance to Langevin phase
NRUNS = 100; % Number of times to run Monte Carlo simulation (>= 1)

% Create Bayes Net from net file
node = makenetwork('net2node.file');

% Number of variables in the Bayes Net
Nvars = length(node);

% Create new ground truth Bayes Net
node = createBN(node, rg);
% Translate BN into sum-product network (circ/meta)
[circ, meta] = cdh_net2circuit(node);
% Number of params is equal to length of meta.varmap
Nparams = length(meta.varmap);
Nnodes = length(circ); % Number of nodes in sum-product network
% Determine locations of free parameters
[fpidx, ~] = cdh_find_free_params(node, meta);
Nfreeparams = length(fpidx);


% For each of the Monte Carlo runs, compute means and variances. On each MC
% run, we generate a new set of observations
mcthetameans = zeros(Nfreeparams, NRUNS);
mcthetavariances = zeros(Nfreeparams, NRUNS);
for mcrun = 1:NRUNS
    tic;
    fprintf(1, '\n---STARTING RUN %d---\n', mcrun);

    % Generate new complete set of Ntrain observations
    val = cdh_sim_observations(node, Ntrain);
    
    % Some of the nodes are unobserved (create mask for which vars are observed
    % over the Ntrain samples).
    obsmask = true(Ntrain, Nvars);
    if PARTIALOBS == true
        obsmask = rand(Ntrain, Nvars) < FRACOBS;
    end

    % Now that we have all observations (including unobserved vars), infer mean
    % and covariance for both complete and incomplete observations
    % RETURNS 0 in cases other than 2-variable Bayes Net case! i.e. doesn't
    % work for those cases, but will not crash
    [mcomp, rcomp] = learnposterior(val);
    [minc, rinc] = learnposterior(val .* obsmask);

    % Apply obsmask to "remove" observations
    val = val .* obsmask;

    % Initialize
    thetadiff = Inf;
    steadystate = 0;
    langevin = false;
    langesamples = 0;
    langesampidx = zeros(1, LANGEMAXSAMP);
    epsilon = 0.001;
    t = 1;
    % Pre-allocate memory to store theta on each iteration
    theta = zeros(Nparams, MAXITER);
    % Initialize first column with random valid probability distribution
    theta(:,t) = rand(Nparams, 1);
    theta(:,t) = cdh_apply_constraints(theta(:,t), meta, node);
    % meta.p is used in the sum-product network inference and must be updated
    meta.p = theta(:,t);

    pevidence = 0;
    grad = zeros(Nparams, 1);

    % MAIN LOOP: perform gradient ascent and Langevin dynamics with
    %   Metropolis-Hastings algorithm to generate samples from the posterior
    %   distribution of theta given the observations
    while (t < MAXITER)
        % epsilon/step size decays unless in Langevin dynamics phase, in which
        % case we hold epsilon constant
        if langevin == false
            epsilon = 0.1/(100 + t);
        end

        % Compute gradient/partial gradient using Nminibatch samples
        [scores_t,grad_t,pe_t] = cdh_bernoulli_fast_gradient(val,fpidx,circ,meta,Nminibatch,t);

        % Gradient ascent step
        th = theta(:,t) + (epsilon/2)*grad_t;

        % Add noise to free params (unless pure gradient ascent enabled)
        if LANGEFLAG == true
            noise = normrnd(0, sqrt(epsilon), Nfreeparams, 1);
            th(fpidx) = th(fpidx) + noise;
        end

        % Threshold free params into interval [0,1]
        % TODO: Generalize to non-binary variables
        for i = 1:Nfreeparams
            if th(fpidx(i)) < 0
                th(fpidx(i)) = 0;
            elseif th(fpidx(i)) > 1
                th(fpidx(i)) = 1;
            end

            th(fpidx(i)+1) = 1 - th(fpidx(i));
        end

        % Determine whether to consider the update part of Langevin dynamics
        % phase
        if LANGEFLAG == true && langevin == false
            % Determine if we are in Langevin dynamics phase, i.e. variance of
            % injected noise greater than variance of (stochastic) gradient step.

            % Compute empirical covariance of scores for free params
            smean = mean(scores_t(fpidx), 2);
            Vs = zeros(Nfreeparams);
            for j = 1:Nminibatch
                Vs = Vs + (scores_t(fpidx,j) - smean)*(scores_t(fpidx,j) - smean)';
            end

            % Compute sample threshold, alpha
            % Note: eigs(Vs, 1) is the max eigenvalue of the empirical
            % covariance matrix
            alpha = epsilon*Ntrain^2/(4*Nminibatch) * eigs(Vs, 1);

            % If alpha << 1, then we are in Langevin dynamics phase
            if alpha < LANGEALPHATHRESH
                if langevin == false
                    fprintf(1, 'langevin f->t, alpha=%4.4f, t=%d\n', alpha, t);
                end
                langevin = true;
            elseif alpha >= 1
                if langevin == true
                    fprintf(1, 'langevin t->f, alpha=%4.4f, t=%d\n', alpha, t);
                end
                langevin = false;
            end

        end

        % Metropolis-Hastings Algorithm to accept or reject the update
        if langevin == true

            % Evaluate multivariate normal at theta{t+1} aka th. Note that the
            % mean of the mvn is at theta{t}+epsilon/2*grad{theta{t}}
            fforw = mvnpdf(th(fpidx)', (theta(fpidx,t) + epsilon/2*grad_t(fpidx))', sqrt(epsilon)*eye(Nfreeparams));

            % Evaluate multivariate normal at theta{t}. The mean of the mvn is
            % at theta{t+1}+epsilon/2*grad{theta{t+1}}. Therefore, we need to
            % compute the gradient with respect to theta{t+1}, aka th.
            % Temporarily set meta.p to th in order to compute new gradient and
            % pevidence for time t+1
            meta.p = th;
            [~,grad_tp1,pe_tp1] = cdh_bernoulli_fast_gradient(val,fpidx,circ,meta,Nminibatch,t);
            fback = mvnpdf(theta(fpidx,t)', (th(fpidx) + epsilon/2*grad_tp1(fpidx))', sqrt(epsilon)*eye(Nfreeparams));     

            % accept probability is (f(th{t}|th{t+1}) *
            % pe_tp1)/(f(th{t+1}|th{t}) * pe_t)
            accept = (fback*pe_tp1)/(fforw*pe_t);

            if accept < 1
                u = rand;
                % Reject update with probability 1 - accept
                if u > accept
                    th = theta(:,t); % Take the previous theta as the new theta
                    fprintf(1, 'MH rejection at time %d\n', t);
                end
            end

            % Update as new sample in Langevin dynamics phase
            langesamples = langesamples + 1;
            langesampidx(langesamples) = t;
            if mod(langesamples, 20) == 0
                fprintf(1, 'Langevin sample %d at time %d\n', langesamples, t);
            end
        end

        % Update meta.p and theta matrix
        if t+1 > size(theta, 2)
            theta = [theta, th];
        else
            theta(:,t+1) = th;
        end

        % th and theta(:,t+1) should now contain the corrected theta for the
        % next iteration
        meta.p = theta(:,t+1);


        if LANGEFLAG == true
            if langesamples >= LANGEMAXSAMP
                break;
            end
        else
            % Compute diff for stopping condition
            thetadiff = norm(theta(:,t+1) - theta(:,t));
            if thetadiff < CONVERGE
                steadystate = steadystate + 1;
                if steadystate >= STOP
                    break;
                end
            end
        end

        t = t + 1;

        if mod(t, 20) == 0
            fprintf(1, '%d\n', t);
        end
    end

    fprintf(1, 'Finished in %d iterations\n', t);

    % Chop off extra columns of theta matrix
    theta = theta(:,1:t);

    % Compute mean and covariance of posterior using Langevin Dynamics
    if LANGEFLAG && langesamples >= LANGEMAXSAMP
        fprintf(1, 'Generated Langevin dynamics estimate of mean/cov\n');
        thetahat = mean(theta(:,langesampidx(LANGEBURNIN+1:end)), 2);
        thetacovhat = cov(theta(:,langesampidx(LANGEBURNIN+1:end))');

        mcthetameans(:,mcrun) = thetahat(fpidx);
        dcov = diag(thetacovhat);
        mcthetavariances(:,mcrun) = dcov(fpidx);
        
        if NRUNS == 1
            figure;
            subplot(1, 2, 1);
            imagesc(thetacovhat);
            title('Langevin Dynamics Covariance Estimate');

            subplot(1, 2, 2);
            [THBAR, THCOV] = cdh_theta_theoretical(node, val, meta);
            imagesc(THCOV);
            colorbar('location','Manual', 'position', [0.93 0.11 0.02 0.81]);
            title('Ground Truth Covariance');
        end
    end 
    
    toc;
end


% Get ground truth probabilities from ground truth BN
pgt = [];
for i = 1:Nvars
    pgt = [pgt; node(i).p(:)];
end
pgt = pgt(fpidx); % Only keep the free parameters
pgt = kron(ones(NRUNS, 1), pgt);

% Convert MC means and variances into subjective opinions
we_so = cdh_gen_subjective_opinions(mcthetameans, mcthetavariances);

% Compute desired confidence bound
[pa,pp] = perfbound_beta_new2(pgt,we_so);
figure;
plot(pa, pp);
title('Desired Confidence Bound Divergence');
xlabel('Desired Confidence');
ylabel('Actual Confidence');

