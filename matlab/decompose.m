function [vo,g] = decompose(node,v)

loc = find(v==1);
ne = length(loc);


vo = zeros(1,size(v,2));
i = 0;
v2 = v;
while sum(v2)>0,
    i = i+1;
    s = find(v2==1,1);
    g{i} = s;
    v2(s) = 0;
    vo(i,s) = 1;
    clean = 1:length(v);
    par = node(s).parents;
    vo(i,par(v(par)==1)) = 1;
    active = par(v(par)==0);
    clean = setdiff(clean,s);
    while ~isempty(active),
       % if i==1, active, end
        pick = active(1);
        g{i} = union(g{i},pick);
        clean = setdiff(clean,pick);
        active = active(2:end);
        par = node(pick).parents;
        vo(i,par(v(par)==1)) = 1;
        active = union(active,par(v(par)==0));
        if v(pick)==0,
            active = union(active,intersect(node(pick).children,clean));
        else
            vo(i,pick) = 1;
            v2(pick) = 0;
        end
    end
end