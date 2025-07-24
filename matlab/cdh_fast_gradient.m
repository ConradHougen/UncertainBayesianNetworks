% Fast gradient computation (reduces number of passes through sum-product
% network)
function [scores,grad,pe] = cdh_fast_gradient(val,fpidx,circ,meta,Nminibatch,t)

    % Nparams is the number of total params (only a fraction are free
    % params)
    Nparams = length(meta.varmap);
    Nnodes = length(circ);
    [Ntrain,Nvars] = size(val);
    
    allgrad = containers.Map; % Initialize cell array for stored gradients
    allpe = containers.Map; % Cell array for stored pevidence
    scores = zeros(Nparams, Nminibatch);
    sidx = 1;
    for i=mod((t-1)*Nminibatch:t*Nminibatch-1, Ntrain) + 1
        v = val(i,:);
        obs = [(1:Nvars)' v'];
        % Unobserved variables are identified by a zero in v. Remove
        % unobserved variable rows in obs
        obsidx = obs(:,2) ~= 0;
        obs = obs(obsidx,:);
        
        % If v is a new observation without a known gradient, run
        % forward/backward pass, otherwise, just use stored gradient
        if ~isKey(allgrad, num2str(v))
                    
            % Forward and backward pass through sum-product network
            [~,m,~,circ] = cdh_inferProb(circ,meta,obs);
            pe = m(Nnodes); % P(evidence) stored in root node of SP-network
        
            % Numerical stability
            if pe < 0.001
                pe = 0.001;
                fprintf(1, 'pevidence=%4.4f, t=%d\n',pe, t);
            end
        
            scores(:,sidx) = m(Nnodes+1:Nnodes+Nparams)/pe;
            
            % Store newly computed gradient and pevidence
            allgrad(num2str(v)) = m(Nnodes+1:Nnodes+Nparams);
            allpe(num2str(v)) = pe;
            
        % Else, we have already computed the gradient and pe for this
        % training example, so just retrieve
        else
            scores(:,sidx) = allgrad(num2str(v))/allpe(num2str(v));
        end
        
        % Increment scores index for next training example
        sidx = sidx + 1;
    end
    
    % Compute gradient of free parameters (located at fpidx)
    % Total derivative = df/dtheta_fp - df/dtheta_NONfp
    grad = sum(scores, 2);
    
    % Prepend 0 to allow ifpidx to be used to divide parameters into sets
    % of parameters which are each a probability distribution 
    ifpidx = [0 setdiff(1:Nparams, fpidx)];
    
    % For each set of parameters, subtract the partial derivative coming
    % from the non-free parameter in the set (i.e. the last parameter in
    % each set)
    for i = 1:length(ifpidx)-1
        grad(ifpidx(i)+1:ifpidx(i+1)-1) = grad(ifpidx(i)+1:ifpidx(i+1)-1) - grad(ifpidx(i+1));
    end

end
