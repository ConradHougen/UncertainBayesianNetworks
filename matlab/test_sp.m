function test_sp(storefile,inferfun)

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
tm0 = zeros(Nruns,1);
cnt = 0;
cnt2 = 0;
for i=1:Nruns,
    waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
    node = nodestore{i};
    obs = getobs(node);
    for j=1:length(node),
        card(j) = node(j).cardinality;
    end
    tic,
    [circ,meta] = net2circuit(node);
    tm0(i) = toc;
    tic,
    [p,v] = feval(inferfun,circ,meta,obs); 
    tm(i) = toc;
    
    for j=1:Nint,
        cnt = cnt+1;
        cd(cnt) = card(interior_loc(j));        
        pp2 = p(interior_loc(j),1:cd(cnt))';
        vv2 = v(interior_loc(j),1:cd(cnt))';
        s = pp2.*(1-pp2)./vv2-1;
        we((cnt2+1):(cnt2+cd(cnt)),:) = [pp2-1./s (1-pp2)-1./s 2./s];
        cnt2 = cnt2+cd(cnt);
    end
end
delete(hd);
we = we(1:cnt2,:);

eval(sprintf('save %s_%s we tm0 tm',inferfun,storefile));
    
    