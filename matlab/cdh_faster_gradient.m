% Fast gradient computation (reduces number of passes through sum-product
% network)
function [scores,grad,pe] = cdh_faster_gradient(val,fpidx,circ,meta,Nminibatch,t)

    ZERO = 1e-10;

    % Nparams is the number of total params (only a fraction are free
    % params)
    Nparams = length(meta.varmap);
    Nnodes = length(circ);
    [Ntrain,Nvars] = size(val);
    
    scores = zeros(Nparams, Ntrain);
    
    [~,m] = inferProbfast(meta,val);
    
    scores = (m(1:Nparams,:)>ZERO).*m(Nnodes+(1:Nparams),:)./(ones(Nparams,1)*m(Nnodes,:));
    pe = exp(sum(log(m(Nnodes,:))));    
    %pe = m(Nnodes,end);
    
    %pe = m(Nnodes,1);
    
    % Iterate over data items to create scores
    %for i = 1:Ntrain
    %    pe = m(Nnodes,i);
    %    scores(:,i) = m(Nnodes+1:Nnodes+Nparams,i)/pe;
    %end
    
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
