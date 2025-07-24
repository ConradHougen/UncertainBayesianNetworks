function [p,v] = we2pv(we)

p = we*[1;0;0.5];

v = p.*(1-p).*we(:,3)./(2+we(:,3));