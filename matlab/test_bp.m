function test_bp(storefile,inferfun,nrounds)

twait = 600; % dump data every 10 minutes
tthresh = twait;

if nargin>2,
    rnd = 1;
else
    rnd = 0;
end

load(storefile);

mxrg = 0; 
for i=1:length(nodestore), 
    for j=1:length(nodestore{i}), 
        mxrg = max(mxrg,nodestore{i}(j).cardinality); 
    end 
end

Nruns = length(nodestore);
node = nodestore{1};
Nint = 0;
interior_loc = [];
for i=1:length(node),
    if isempty(node(i).value);
        Nint = Nint+1;
        interior_loc = [interior_loc i];
    end
end

hd = waitbar(0);

we = zeros(Nruns*Nint*mxrg,3);
tm = zeros(Nruns,1);
cnt = 0;
cnt2 = 0;

for i=1:Nruns,
    waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
    node = nodestore{i};
    
    for j=1:length(node),
        card(j) = node(j).cardinality;
    end
    if rnd,
        tic,
        node = feval(inferfun,node,nrounds); 
        tm(i) = toc;
        if sum(tm)>tthresh,
            tthersh = tthersh + twait;
            eval(sprintf('save %s_%d_%s we tm i',inferfun,nrounds,storefile));
        end
    else
        tic,
        node = feval(inferfun,node); 
        tm(i) = toc;
    end
    
    for j=1:Nint,
        cnt = cnt+1;
        cd(cnt) = card(interior_loc(j));        
        pp = node(interior_loc(j)).me';
        vv = diag(node(interior_loc(j)).mec);
        s = pp.*(1-pp)./vv-1;
        we((cnt2+1):(cnt2+cd(cnt)),:) = [pp-1./s (1-pp)-1./s 2./s];
        cnt2 = cnt2+cd(cnt);
    end
end
delete(hd);
we = we(1:cnt2,:);

if rnd,
    eval(sprintf('save %s_%d_%s we tm',inferfun,nrounds,storefile));
else
    eval(sprintf('save %s_%s we tm',inferfun,storefile));
end
    
    