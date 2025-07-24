function [m,cov,theta0] = mh_fix(meta,val)

M = 1000;

Nparam = length(meta.varmap);

map = meta.valuemap(:,1:Nparam);
card = max(map,[],2);

N = size(map,1); 
Ntrain = size(val,1);


% nsee = sum(val>0,2);
% valc = val(nsee==N,:);
% 
% cs = [1; cumprod(card(1:end-1))];
% 
% theta0 = ones(Nparam,1);
% for i=1:Nparam,
%     idx = max(meta.valuemap(:,i)-1,0)'*cs;
%     vc = valc;
%     vc(:,meta.valuemap(:,i)==0) = 0;
%     theta0(i) = theta0(i)+sum(max(vc-1,0)*cs==idx);
% end

v = reshape(meta.valuemap(:,1:Nparam),[N 1 Nparam]); 
v = repmat(v,[1 Ntrain 1]);

na = sum(meta.valuemap(:,1:Nparam)>0,1);
na = reshape(na,[1 1 Nparam]);
na = repmat(na,[1 Ntrain 1]);

valn = val;
valn(val==0) = NaN;
hits = sum(repmat(valn',[1 1 Nparam])==v,1)==na;
hits = double(reshape(hits,[Ntrain Nparam]));

theta0 = sum(hits,1)'+1;



    

cnt = 0;
mask = zeros(Nparam,Nparam);
for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        cnt = cnt+card(i);
    end
end

r = gamrnd(theta0*ones(1,M),1); 
r = r./(mask*r);

g = exp(sum((theta0*ones(1,M)-1).*log(r),1));
post = zeros(1,M); 

for i=1:M,
    meta.p = r(:,i);
    post(i) = inferProbfor(meta,val);
end

% A = min(post(2:M).*g(1:M-1)./post(1:M-1)./g(2:M),1);
% 
% a = rand(1,M-1)<A;

for i=1:M-1,
    a = rand(1,1)<min(post(i+1)*g(i)/post(i)/g(i+1),1); 
    if a==0,
        r(:,i+1) = r(:,i);
        post(i+1) = post(i);
        g(i+1) = g(i);
    end
end


%r = r(:,1000:end);
%M = size(r,2);
m = mean(r,2);
r = r-m*ones(1,M);
cov = (r*r')/M;

