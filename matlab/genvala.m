function vala = genvala(meta)

Nparam = length(meta.varmap);

map = meta.valuemap(:,1:Nparam);
card = max(map,[],2)+1;

N = size(map,1);

v = zeros(1,N);

Nv = prod(card);

vala = zeros(Nv,N);
for i=2:Nv,
    id = N;
    while (id>0),
        v(id) = v(id)+1;
        if v(id)==card(id),
            v(id) = 0;
            id = id-1;
        else
            id = 0;
        end
    end
    vala(i,:) = v;
end
