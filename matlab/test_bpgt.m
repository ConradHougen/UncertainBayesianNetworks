function test_bpgt(storefile,inferfun)

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

pgt = zeros(Nruns*Nint*mxrg,1);
tm = zeros(Nruns,1);
cnt = 0;
cnt2 = 0;
for i=1:Nruns,
    waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
    node = nodestore{i};
    
    for j=1:length(node),
        card(j) = node(j).cardinality;
    end
    tic,
    node = feval(inferfun,node); 
    tm(i) = toc;
    
    for j=1:Nint,
        cnt = cnt+1;
        cd(cnt) = card(interior_loc(j));        
        pgt((cnt2+1):(cnt2+cd(cnt))) = node(interior_loc(j)).pe';
        cnt2 = cnt2+cd(cnt);
    end
end
delete(hd);
pgt = pgt(1:cnt2);

eval(sprintf('save %s_%s pgt tm',inferfun,storefile));
    
    