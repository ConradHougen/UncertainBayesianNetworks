node = makenetwork('net2node.file');

Ntrain = 100;
Ncomp = 20;
Nruns = 200;


pgt = [];
we = [];
wei = [];
wm = [];

work.val = zeros(100,2,Nruns);
work.cov = zeros(6,6,Nruns);
work.pest = zeros(6,Nruns);
work.pm = zeros(6,Nruns);
work.covm = zeros(6,6,Nruns);
work.pgt = zeros(6,Nruns);


for i=1:Nruns,
    
    disp(num2str(i));
    
    node = createBN(node,[2 2]);
    [node,val] = createSLBN(node,Ntrain);
    [circ,meta] = net2circuit(node);
    
    work.node{i} = node;
    
    val = val+1;
    
    work.val(:,:,i) = val;
    
    val((Ncomp+1):end,1) = 0;
    

    [pem,cov,mask] = em_fast(circ,meta,val);
    
    work.cov(:,:,i) = cov;
    work.pest(:,i) = pem;
    work.pgt(:,i) = meta.p;
    
    [pm,covm] = mm(circ,meta,val);
    
    work.covm(:,:,i) = covm;
    work.pm(:,i) = pm;
    
    
    meta.m = pem;
    meta.cov = cov;
    obs = [2 val(1,:)];
    pg = inferProb(circ,meta,obs); 
    [p,v] = inferSOProb(circ,meta,obs);
    pg = pg(1,:);
    pgt = [pgt; pg(:)];
    p = p(1,:); 
    v = v(1,:);
    p = p(:);
    v = v(:);
    s = p.*(1-p)./v-1;
    s = max([1./p 1./(1-p) s],[],2);
      
%     if sum(s<0)>0,
%         break
%     end
    we = [we; [p-1./s (1-p)-1./s 2./s]];
    
    meta.cov = cov.*mask;
    [p,v] = inferSOProb(circ,meta,obs);
    p = p(1,:); 
    v = v(1,:);
    p = p(:);
    v = v(:);
    s = p.*(1-p)./v-1;
    s = max([1./p 1./(1-p) s],[],2);
    
    wei = [wei; [p-1./s (1-p)-1./s 2./s]];
    
    meta.m = pm;
    meta.cov = covm;
    [p,v] = inferSOProb(circ,meta,obs);
    p = p(1,:); 
    v = v(1,:);
    p = p(:);
    v = v(:);
    s = p.*(1-p)./v-1;
    s = max([1./p 1./(1-p) s],[],2);
    wm = [wm; [p-1./s (1-p)-1./s 2./s]];
    
end

[pa,pp] = perfbound_beta_new2(pgt,we);
[pa2,pp] = perfbound_beta_new2(pgt,wm);
hd = plot(pp,pa,'b-',pp,pa2,'r-',pp,pp,'k:');

%save work_em.mat



    