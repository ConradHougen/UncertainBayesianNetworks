function [nodeout,out] = inferMonteBP(node,nrounds)

if nargin<2,
    nrounds = 100;
end

N = length(node);

out = struct('alpha',cell(1,N),'pe',cell(1,N));
for j=1:N,
    w = node(j).w;
    l = size(w,2)-1;
    out(j).alpha = l*w(:,1:end-1)./(w(:,end)*ones(1,l))+1;
    out(j).pe = zeros(nrounds,l);
end

nodeout = node;

for i=1:nrounds,
    for j=1:N,
        p = gamrnd(out(j).alpha,1);
        p = p./(sum(p,2)*ones(1,node(j).cardinality));
        node(j).p = p';
    end 
    node = prob_infer(node);
    for j=1:N,
        out(j).pe(i,:) = node(j).pe;
    end
end

for i=1:N,
    if ~isempty(nodeout(i).value),
        nodeout(i).we = zeros(1,nodeout(i).cardinality);
        nodeout(i).we(nodeout(i).value+1) = 1;
        nodeout(i).wb = zeros(nodeout(i).cardinality,3);
        nodeout(i).wb(:,2) = 1;
        nodeout(i).wb(nodeout(i).value+1,:) = [1 0 0];
    else
        m = mean(out(i).pe,1);
        v = (out(i).pe'*out(i).pe)/nrounds;
%         v = mean(out(i).pe.^2,1);
%         s = sum(m.*(1-m).*(m-v))/sum(m.*(1-m).*(v-m.^2));
% %         s = (m-v)./(v-m.^2);
% %         s = s(2);
%         s = max([s 1./m]);
%         nodeout(i).we = [m-1/s nodeout(i).cardinality/s];
%         s = (m-v)./(v-m.^2);
%         s = 1./s';
%         nodeout(i).wb = [m'-s 1-m'-s 2*s];
        nodeout(i).me = m;
        nodeout(i).mec = v - m'*m;
    end
end
    
    