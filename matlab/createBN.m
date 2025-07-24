function node = createBN(node,rg)

if nargin<2,
    rg = [2 5];
end

% This function creates a random multinomial BN for a given DAG

N = length(node);

cardinality = randi(rg, 1, N); % creates a row length N with the cardinalities of the nodes
                                  % as randomly generated numbers between 2 and 10
                                  % and then if |X|=k
                                  % we will assume the values it takes
                                  % are 0,1,...,k-1

for i=1:N,
    node(i).cardinality = cardinality(i);  %Added by LMK
    node(i).pcard = cardinality(node(i).parents); %Added by LMK
    nparents = length(node(i).parents); % number of parents of the node i
    card_p = prod(cardinality(node(i).parents));% product of cardinalities of 
                                                         % the parents of i
    node(i).p = []; % initiate the conditional probability matrix of the node i
    for m=1:card_p,
        A=rand(cardinality(i),1); %a column of random numbers size cardinality(i)
        node(i).p = [node(i).p, A/sum(A)]; 
        % node(i).p contains one column of random multinomial distribution
        % of cardinality(i) for each setting of the parents of the node i 
    end
                                    
end

%%%LMK - changed variable cardinality(i).p to card_p
%%% Added lines 14-15 to store cardinality of each variable node