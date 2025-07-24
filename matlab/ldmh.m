function [phat,chat,theta] = ldmh(out,meta,val,langeparams,i)

Nparams = length(meta.varmap);
MAXITER = langeparams.MAXITER;
LANGEMAXSAMP = langeparams.LANGEMAXSAMP;
LANGEBURNIN = langeparams.LANGEBURNIN;
LANGEALPHATHRESH = langeparams.LANGEALPHATHRESH;
Ntrain = langeparams.Ntrain;
Nminibatch = langeparams.Nminibatch;
fpidx = langeparams.fpidx;
ifpidx = langeparams.ifpidx;
circ = langeparams.circ;
Nfreeparams = length(fpidx);

langevin = false; % flag to indicate entering langevin dynamics
langesamples = 0;
langesampidx = zeros(1, LANGEMAXSAMP);
epsilon = 0.001; % Initial step size/variance
t = 1;

theta = zeros(Nparams, MAXITER);
theta(:,t) = rand(Nparams, 1);
theta(:,t) = cdh_apply_constraints(theta(:,t),meta,out.node(i,:));

metalocal = meta;
metalocal.p = theta(:,t);

% LD loop to converge to parameter and covariance estimates
while t < MAXITER
    % Epsilon decay until Langevin dynamics phase
    if langevin == false
        epsilon = 0.1/(100 + t);
    end

    % Compute gradient/partial gradient using Nminibatch samples
    [scores_t,grad_t,pe_t] = cdh_faster_gradient(val,fpidx,circ,metalocal,Nminibatch,t);

    % Gradient ascent step 
    grad_t(ifpidx(2:end)) = 0; % (only free params should be updated)
    th = theta(:,t) + (epsilon/2)*grad_t;

    % Add noise to free params
    noise = normrnd(0, sqrt(epsilon), Nfreeparams, 1);
    th(fpidx) = th(fpidx) + noise;

    % Apply KKT conditions
    th = cdh_kkt_constrain(th, ifpidx);

    % Determine whether to consider the update part of Langevin dynamics
    % phase
    if langevin == false
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
        metalocal.p = th;
        [~,grad_tp1,pe_tp1] = cdh_faster_gradient(val,fpidx,circ,metalocal,Nminibatch,t);
        fback = mvnpdf(theta(fpidx,t)', (th(fpidx) + epsilon/2*grad_tp1(fpidx))', sqrt(epsilon)*eye(Nfreeparams));     

        % accept probability is (f(th{t}|th{t+1}) *
        % pe_tp1)/(f(th{t+1}|th{t}) * pe_t)
        accept = (fback*pe_tp1)/(fforw*pe_t);

        if accept < 1
            u = rand;
            % Reject update with probability 1 - accept
            if u > accept
                th = theta(:,t); % Take the previous theta as the new theta
                %fprintf(1, 'MH rejection at time %d\n', t);
            end
        end

        % Update as new sample in Langevin dynamics phase
        langesamples = langesamples + 1;
        langesampidx(langesamples) = t;
        if mod(langesamples, 50) == 0
            fprintf(1, 'Langevin sample %d at time %d\n', langesamples, t);
        end
    end

    % Update theta matrix
    if t+1 > size(theta, 2)
        theta = [theta, th];
    else
        theta(:,t+1) = th;
    end

    % th and theta(:,t+1) should now contain the corrected theta for the
    % next iteration
    metalocal.p = theta(:,t+1);

    if langesamples >= LANGEMAXSAMP
        break;
    end


    t = t + 1;

    if mod(t, 20) == 0
        fprintf(1, '%d\n', t);
    end
end

fprintf(1, 'Finished in %d iterations\n', t);

% Chop off extra columns of theta matrix
theta = theta(:,1:t);

phat = mean(theta(:,langesampidx(LANGEBURNIN+1:end)),2);
chat = cov(theta(:,langesampidx(LANGEBURNIN+1:end))');

end