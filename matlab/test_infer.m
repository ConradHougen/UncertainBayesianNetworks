    
% This is the main file that calls code to perform regular belief 
% propagation over the ground truth conditional probabilities and 
% subjective belief propagation over the conditional opinions.


% Ntrain: Number of full observations of instantiations of variable 
% values to form the conditional opinions in the subjective BN.

% Nnetwork: Number of ground truth Bayesian networks with the given 
% structure.  A specific network has specific conditional 
% probability values.

% Nmonte: Number of different subjective Bayesian networks 
% formed from each ground truth Bayesian network.
clear;
useBP = 0;  %Calculate BP and Second order Extensions of BP
useLBP = 1; %Calulate Loopy BP and second order extensions of Loopy BP
     
rg = [2 2];   %Cardinality range for the varibles
mxrg = rg(2);

Ntrain = 50;  %Number of instantiations for training
Nmonte = 100; %100; %Number of SLBNs created per network
Nnetworks = 10; %10; %Number of random Bayesian netowrks

%load state.mat
%rand('state',st);

% For each ground truth BN, 10 SBNs are formed 


Nruns = Nmonte*Nnetworks;

hd = waitbar(0);


node = makenetwork('27nodenet.file');

% the lists of active and non-active (interior) nodes and their lengths:

N = length(node);
end_loc = getends(node);%the starting list of active nodes - those that have one parent or
%child

%end_loc = [1 5];

Nends = length(end_loc);
interior_loc = setdiff((1:N)',end_loc);% the list of non-active nodes
Nint = length(interior_loc);

% pgt = zeros(Nruns*Nint,mxrg);
% we = zeros(Nruns*Nint,mxrg+1);
pgt = zeros(Nruns*Nint*mxrg,1);
pgt2 = zeros(Nruns*Nint*mxrg,1);
pgt3 = zeros(Nruns*Nint*mxrg,1);
pgtLBP = zeros(Nruns*Nint*mxrg,1);
we_so = zeros(Nruns*Nint*mxrg,3);
we_mnvar = zeros(Nruns*Nint*mxrg,3);
we_sobp = zeros(Nruns*Nint*mxrg,3);
we_sobp2 = zeros(Nruns*Nint*mxrg,3);
we_solbp = zeros(Nruns*Nint*mxrg,3);

tm0 = zeros(Nruns,1);
tm_so = zeros(Nruns,1);
tm_fo = zeros(Nruns,1);
tm_fobp= zeros(Nruns,1);
tm_fobp2= zeros(Nruns,1);
tm_folbp = zeros(Nruns,1);
tm_sobp = zeros(Nruns,1);
tm_sobp2 = zeros(Nruns,1);
tm_solbp = zeros(Nruns,1);
tm_slbn = zeros(Nruns,1);
tm_mnvar = zeros(Nruns,1);
% wem = we;
%wec = we;
cd = zeros(Nruns*Nint,1);

card = zeros(1,N);

obs = zeros(Nends,2);

cnt = 0;
cnt2 = 0;
for i=1:Nruns,
    %waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
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
    
    
    
    tic,
    [circ,meta] = net2circuit(node);
    tm0(i) = toc;
    tic,
    p1 = inferProb(circ,meta,obs);
    tm_fo(i) = toc;
    
   
    tic,
    [p_so,v_so] = inferSOProb(circ,meta,obs);
    tm_so(i) = toc;
    tic,
    %[p_mnvar,v_mnvar] = meanvar(circ,meta,obs);
    tm_mnvar(i) = toc;
    
    if useBP==1,
        tic,
        node2 = prob_infer(node);
        tm_fobp(i) = toc;
        tic,
        node2 = inferSO(node2);
        tm_sobp(i) = toc;
        tic,
        node3 = prob_infer4(node);
        tm_fobp2(i) = toc;
        tic,
        node3 = inferSO4(node3);
        tm_sobp2(i) = toc;
    end
    
    % Loopy BP
    if useLBP == 1
        tic,
        nodeLBP_pi = prob_infer_loopy(node,obs);
        tm_folbp(i) = toc;
        tic,
        nodeLBP = inferSOloopy(node,obs);
        tm_solbp(i) = toc;
        %tic;
        %nodeLBP = inferSO(nodeLBP);
        %tm_solbp(i) = toc;
    end
    

    for j=1:Nint,
        cnt = cnt+1;
        cd(cnt) = card(interior_loc(j));        
        pgt((cnt2+1):(cnt2+cd(cnt))) = p1(interior_loc(j),1:cd(cnt));
        pp2 = p_so(interior_loc(j),1:cd(cnt))';
        vv2 = v_so(interior_loc(j),1:cd(cnt))';
        s = pp2.*(1-pp2)./vv2-1;
        we_so((cnt2+1):(cnt2+cd(cnt)),:) = [pp2-1./s (1-pp2)-1./s 2./s];
        
%         pp2 = p_mnvar(interior_loc(j),1:cd(cnt))';
%         vv2 = v_mnvar(interior_loc(j),1:cd(cnt))';
%         s = pp2.*(1-pp2)./vv2-1;
%         we_mnvar((cnt2+1):(cnt2+cd(cnt)),:) = [pp2-1./s (1-pp2)-1./s 2./s];
%         
        if useBP==1,
            pgt2((cnt2+1):(cnt2+cd(cnt))) = node2(interior_loc(j)).pe';
            pp2bp = node2(interior_loc(j)).me';
            vv2 = diag(node2(interior_loc(j)).mec);
            s = pp2bp.*(1-pp2bp)./vv2-1;
            we_sobp((cnt2+1):(cnt2+cd(cnt)),:) = [pp2bp-1./s (1-pp2bp)-1./s 2./s];
            pgt3((cnt2+1):(cnt2+cd(cnt))) = node3(interior_loc(j)).pe';
            pp2bp = node3(interior_loc(j)).me';
            vv2 = diag(node3(interior_loc(j)).mec);
            s = pp2bp.*(1-pp2bp)./vv2-1;
            we_sobp2((cnt2+1):(cnt2+cd(cnt)),:) = [pp2bp-1./s (1-pp2bp)-1./s 2./s];   
        end
        
%         if useLBP==1
%             pgtLBP((cnt2+1):(cnt2+cd(cnt))) = nodeLBP(interior_loc(j)).pe';
%             pp2lbp = nodeLBP(interior_loc(j)).me';
%             vv2lbp = diag(nodeLBP(interior_loc(j)).mec);
%             s = pp2lbp.*(1-pp2lbp)./vv2lbp-1;
%             we_solbp((cnt2+1):(cnt2+cd(cnt)),:) = [pp2lbp-1./s (1-pp2lbp)-1./s 2./s];
%             pgtLBP((cnt2+1):(cnt2+cd(cnt))) = nodeLBP(interior_loc(j)).pe';
%             pp2lbp = nodeLBP(interior_loc(j)).me';
%             vv2lbp = diag(nodeLBP(interior_loc(j)).mec);
%             s = pp2lbp.*(1-pp2lbp)./vv2lbp-1;
%             we_solbp((cnt2+1):(cnt2+cd(cnt)),:) = [pp2lbp-1./s (1-pp2lbp)-1./s 2./s];
%         end
        
               
        cnt2 = cnt2+cd(cnt);        
%         wem(cnt,1:card(interior_loc(j))+1) = nodem(interior_loc(j)).we;
    %    wec(cnt,:) = nodec(interior_loc(j)).we;
        %if isnan(wec(cnt,1)),
     %   if sum(wec(cnt,:)<0)>0,
      %      break,
       % end
    end
%     if sum(wec(cnt,:)<0)>0,
%         break
%     end
end
delete(hd);

pgt = pgt(1:cnt2);
pgt2 = pgt2(1:cnt2);
pgt3 = pgt3(1:cnt2);
pgtLBP = pgtLBP(1:cnt2);
we_so = we_so(1:cnt2,:);
we_mnvar = we_mnvar(1:cnt2,:);
we_sobp = we_sobp(1:cnt2,:);
we_sobp2 = we_sobp2(1:cnt2,:);
we_solbp = we_solbp(1:cnt2,:);


%[we3,pgt3] = sortoutresults(we,pgt,cd);

[pa,pp] = perfbound_beta_new2(pgt,we_so);
pa2 = perfbound_beta_new2(pgt,we_mnvar);

% [wem2,pgt2] = sortoutresults(wem,pgt,cd);
% 
%pam = perfbound_beta_new2(pgt3,we3);

% [pac,pp] = perfbound_beta(pgt,wec);
% 
% plot(pp,pa,'b-',pp,pac,'r-',pp,pp,'k:');

plot(pp,pa,'b-',pp,pa2,'r-',pp,pp,'k:');

%plot(pp,pa,'b-',pp,pp,'k:');
    
nmax = @(p1,p2) max(abs(p1-p2)./max(p1,p2));

[p1,v1] = we2pv(we_so);
[p2,v2] = we2pv(we_mnvar);
[p3,v3] = we2pv(we_sobp);
[p4,v4] = we2pv(we_solbp);

ep12 = nmax(p1,p2);  %Calclaulate max relative errors for mean and variance estimates
ep13 = nmax(p1,p3);
ep23 = nmax(p2,p3);

ev12 = nmax(v1,v2);
ev13 = nmax(v1,v3);
ev23 = nmax(v2,v3);
    
