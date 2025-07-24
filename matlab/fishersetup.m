function metafisher = fishersetup(meta,node)

N = size(meta.valuemap,1);

Nparam = length(meta.varmap);
map = meta.valuemap(:,1:Nparam);
card = max(map,[],2);

cnt = 0;
cnt2 = 0;
mask = zeros(Nparam,Nparam);
for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        mask2(cnt2+1:cnt2+card(i)-1,cnt+1:cnt+card(i)) = [eye(card(i)-1) -ones(card(i)-1,1)];
        cnt = cnt+card(i);
        cnt2 = cnt2+card(i)-1;
    end
end

vala = genvala(meta);

%idvala = (vala>0)*(2.^(N-1:-1:0)');

Nc = 2^N-1;

for i=1:N,
    loc = find(meta.varmap==i);
    loc2 = loc(1);
    v = zeros(Nparam,1);
    v(loc) = 1;
    v = (mask2>0)*v;
    loc = find(v==1);
    v = zeros(1,N);
    v(meta.valuemap(:,loc2)>0) = 1;
    id = v*(2.^(N-1:-1:0)');
    jcomputed(id) = 1;
    metafisher.loc{id} = loc; 
end




for i = 1:Nc,
    if i>0,
        v = dec2bin(i,N)-48;
        [vo,g] = decompose(node,v);
        id = vo*(2.^(N-1:-1:0)');
        metafisher.dec{i} = id;
        for j=1:size(vo,1),
            if ~jcomputed(id(j)),
                 loc = [];
                 for k=1:length(g{j}),
                     loc = union(loc,find(meta.varmap==g{j}(k)));
                 end
                 metafisher.loc{id(j)} = loc;  
            end
        end
    end
end

metafisher.vala = vala;