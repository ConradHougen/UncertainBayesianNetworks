function [node, val] = createSLBN(node,ntrain)

%for a given BN (node) and given number of observations (ntrain), this
%function simulates ntrain complete observations of the nodes and creates SLBN 
%based on them


N = length(node); %the size of the network
queue = 1:N; %the initial row of non-activated nodes
activated = zeros(N,1); %no node is active at the beginning
order = []; % the row of nodes in the order of activation


while ~isempty(queue),
    for i=queue, % for i=1 to N
        if sum(activated(node(i).parents))==length(node(i).parents),% if all the parents of a node are active
            activated(i) = 1;  % then activate the node
            order = [order i]; 
        end
    end
    queue = setdiff(queue,order);
end

val = zeros(ntrain,N); %the initial matrix of observations, one column per each node

%generating observations at each node based on the observations of its parents and the 
%conditional ditributions at the node in the BN (stohastic simulation?): 

%node(i).p is a matrix determining the conditional distributions at i,
%one column per each setting of the parents
%node(i).parents is the row of parents of i in an increasing order


for i = order,
    np = length(node(i).parents); %number of parents of the node i
    if np==0, % node(i).p is just a column in this case
       
            val(:,i) = discretesample(node(i).p,ntrain);%random samples of X_i go in the ith column of val
        
    else
            A = node(i).pcard; %A = cardinality(node(i).parents);
            B=[A(2:length(A)), 1]; 
            %for ii = 1:ntrain,
                val(:,i) = discretesample(node(i).p(:, val(:,node(i).parents)*cumprod(B,'reverse')'+1));
            %end
            % the k-th sample for node i determined according to its parent
            % samples: the corresponding column from node(i).p is chosen
            % according to the entries in val for the
            % parents for sampling with discretesample 
                          
                                 
    end
end


% counting the observations in val and calculating the corresponding opinions:

for i=1:N,
    np = length(node(i).parents);
    count = zeros(1,node(i).cardinality);
    if np==0,
        for k=0:(node(i).cardinality-1),
            count(k+1) = sum(val(:,i) == k); %number of appearances of k in the val column
            %corresponding to node i 
        end
        node(i).w = [count  node(i).cardinality]/(sum(count)+node(i).cardinality); %opinion at the node i 
    else
        A = node(i).pcard;
        B = [A(2:length(A)) 1];
        idx = val(:,node(i).parents)*cumprod(B,'reverse')'+1;
        %idx-column that determines which element of node(i).p should be
        %taken at each of the rounds
        node(i).w = zeros(prod(node(i).pcard),node(i).cardinality+1);
        for j=1:prod(node(i).pcard), %cardinality(node(i).parents)),
            for k=0:(node(i).cardinality-1),
            count(k+1) = sum(val(idx==j,i)== k);%in the i-th column, counts the number of rows where idx=j
            end
            node(i).w(j,:) = [count node(i).cardinality]/(sum(count)+node(i).cardinality);%we take W=1/|X| 
            %end
            %the opinion matrix has one row per each setting of the parents
        end
    end
end


            
     