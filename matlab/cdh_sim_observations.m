% for a given BN (node), simulate ntrain complete observations, then
% reduces the number of returned observations based on the frac variable.
% E.g., ntrain = 100 and frac = 0.5 means that we will simulate 100
% observations, but then we will assume that (1 - 0.5) * 100 = 50 of those
% values are not actually observed.

function val = cdh_sim_observations(node, nsamples)

N = length(node); % the size of the network
queue = 1:N; % the initial row of non-activated nodes
activated = zeros(N,1); % no node is active at the beginning
order = []; % the row of nodes in the order of activation

while ~isempty(queue)
    for i = queue
        % If all the parents of a node are active, then activate the
        % node
        if sum(activated(node(i).parents)) == length(node(i).parents)
            activated(i) = 1;
            order = [order i];
        end
    end
    queue = setdiff(queue, order);
end

% Initial matrix of observations (N columns for N nodes, ntrain rows
% for ntrain observations per node)
val = zeros(nsamples, N);


%generating observations at each node based on the observations of its parents and the 
%conditional ditributions at the node in the BN (stohastic simulation?): 

%node(i).p is a matrix determining the conditional distributions at i,
%one column per each setting of the parents
%node(i).parents is the row of parents of i in an increasing order

for i = order
    np = length(node(i).parents); % np = number of parents of the node i
    % If the node has no parents, node(i).p has only one column, based
    % on the values that node i can take, only
    if np == 0
        %random samples of X_i go in the ith column of val
        val(:,i) = discretesample(node(i).p, nsamples); 

    else
        A = node(i).pcard;
        B = [A(2:end), 1];
        val(:,i) = discretesample(node(i).p(:, val(:,node(i).parents)*cumprod(B,'reverse')'+1));            
    end

end
    
% +1 to change range to [1,card] instead of [0,card-1]
val = val + 1;
    
    
    
end