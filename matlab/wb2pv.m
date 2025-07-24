function [p,v] = wb2pv(node)

if 

n = length(node);
c = zeros(n,1);
for i=1:n,
    c(i) = node(i).cardinality;
end

p = NaN*ones(n,max(c));
v = p;

for i=1:n,
    m = node(i).wb*[1;0;0.5];
    v(i,1:c(i)) = (m.*(1-m).*node(i).wb(:,3)./(node(i).wb(:,3)+c(i)))';
    p(i,1:c(i)) = m';
end
    
