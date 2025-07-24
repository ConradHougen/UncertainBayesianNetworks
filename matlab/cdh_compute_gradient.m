% Computes full or partial gradient, given observations in val matrix,
% sum-product network information stored in circ and meta, and the number
% of observations to use in this mini batch for the full or partial
% gradient calculation
%
% Returns: a score matrix which contains in each column the "gradient"
% computed using a single sample, the final partial or full gradient in the
% grad variable, and the joint probability of the given observations
% (probability of the evidence) in pe.
function [scores,grad,pe] = cdh_compute_gradient(val,fpidx,circ,meta,Nminibatch,t)

    % Nparams is the number of total params (only a fraction are free
    % params)
    Nparams = length(meta.varmap);
    Nnodes = length(circ);
    [Ntrain,Nvars] = size(val);
    
    grad = zeros(Nparams, 1);
    scores = zeros(Nparams, Nminibatch);
    sidx = 1;
    for i=mod((t-1)*Nminibatch:t*Nminibatch-1, Ntrain) + 1
        v = val(i,:);
        obs = [(1:Nvars)' v'];
        
        % Unobserved variables are identified by a zero in v. Remove
        % unobserved variable rows in obs
        obsidx = obs(:,2) ~= 0;
        obs = obs(obsidx,:);
        
        % Forward and backward pass through sum-product network
        [~,m,~,circ] = cdh_inferProb(circ,meta,obs);
        pe = m(Nnodes); % P(evidence) stored in root node of SP-network
        
        % Numerical stability
        if pe < 0.001
            pe = 0.001;
            fprintf(1, 'pevidence=%4.4f, t=%d\n',pe, t);
        end
        
        scores(:,sidx) = m(Nnodes+1:Nnodes+Nparams)/pe;
        
        % Compute gradient of free parameters (located at fpidx)
        % Total derivative = df/dtheta_xk - df/dtheta_NOTxk
        % TODO: Generalize to non-binary variables
        %grad(fpidx) = grad(fpidx) + (scores(fpidx,sidx) - scores(fpidx+1,sidx));
        
        % Increment scores index for next training example
        sidx = sidx + 1;
    end
    
    grad(fpidx) = sum(scores(fpidx,:), 2) - sum(scores(fpidx+1,:), 2);

end
