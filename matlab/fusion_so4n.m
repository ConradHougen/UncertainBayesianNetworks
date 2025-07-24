function [mo,covo] = fusion_so(mi,covi)

n = length(mi);
ne = length(mi{1});
m = zeros(n,ne);
for i=1:n,
    m(i,:) = mi{i};
end
mo = prod(m,1);
mo = mo/sum(mo);

covo = zeros(ne,ne);
for i=1:n,
    H = mo./mi{i};
    loc = find(isnan(H));
    if ~isempty(loc),
        H(loc) = prod(m([1:i-1 i+1:n],loc));
    end
    H = diag(H)-mo'*H;
    %H = (diag(mo)-mo'*mo)./(ones(ne,1)*mi{i});
    covo = covo + H*covi{i}*H';
end

