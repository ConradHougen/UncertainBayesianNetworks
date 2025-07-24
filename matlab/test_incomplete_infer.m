% Main file

% Given a Bayesian network and incomplete training data, we would like to
% estimate the BN parameters (means and covariances).

% Method:
% 1. Input: Bayes Net Structure in a net file
% 2. Generate a ground truth BN (specific probability assignments)
% 3. Given the ground truth BN, simulate incomplete observations of the
%       nodes
% 4. Generate the sum product network associated with the BN
% 5. Forward/backward pass through sum product net to compute partial
%       derivatives (aka conditional probabilities)
% 6. Compute gradient using partial derivatives, then use Langevin Dynamics
%       to infer means/variances of parameters



rg = [3 3];   %Cardinality range for the variables

Nnetworks = 10; % Number of random Bayesian networks

% 1. Read in BN structure from net file
node = makenetwork('net1.file');

N = length(node);
end_loc = getends(node);

Nends = length(end_loc);
interior_loc = setdiff((1:N)',end_loc); % list of non-active nodes
Nint = length(interior_loc);


card = zeros(1, N);

% generate a new ground truth BN
node = createBN(node, rg);
for j=1:N
    card(j) = node(j).cardinality;
end
    
loc = getends(node);
   
% Generate sum-product network from ground truth BN
[node,val] = createSLBN(node,100);
[circ, meta] = net2circuit(node);

obs = [(1:9)' v'];
[p, m] = inferProb(circ, meta, obs);
    
    

    