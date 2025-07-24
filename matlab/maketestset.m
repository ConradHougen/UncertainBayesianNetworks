rg = [2 4];   %Cardinality range for the varibles
mxrg = rg(2);

Ntrain = 10;  %Number of instantiations for training
Nmonte = 100; %100; %Number of SLBNs created per network
Nnetworks = 10; %10; %Number of random Bayesian netowrks

load state.mat
rand('state',st);

% For each ground truth BN, 10 SBNs are formed 


Nruns = Nmonte*Nnetworks;

hd = waitbar(0);


node = makenetwork('net3.file');

% the lists of active and non-active (interior) nodes and their lengths:

N = length(node);
end_loc = getends(node);%the starting list of active nodes - those that have one parent or
%child

%end_loc = [1 5];

Nends = length(end_loc);
interior_loc = setdiff((1:N)',end_loc);% the list of non-active nodes
Nint = length(interior_loc);

for i=1:Nruns,
    waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
    if mod(i-1,Nmonte)==0, 
        node = createBN(node,rg);% a new ground truth BN is created every 10th run
        for j=1:N,
            card(j) = node(j).cardinality;
        end
    end
    loc = getends(node);
    [node,val] = createSLBN(node,Ntrain);

    for j=1:Nends,
        value = sum(rand(1)>(1:card(end_loc(j))-1)/card(end_loc(j)));
        obs(j,:) = [end_loc(j) 1+value]; %random observation for the end nodes
        node(end_loc(j)).value = value;
    end
    
    nodestore{i} = node;
    
end

delete(hd);


    
    




