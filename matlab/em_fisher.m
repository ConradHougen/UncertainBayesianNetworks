function [pout,cov,mask,mask2,mask3] = em_fisher(meta,val,vala,p0)

diagnostic = 0;

Nparam = length(meta.varmap);

map = meta.valuemap(:,1:Nparam);
card = max(map,[],2);

N = size(map,1); 
Nnodes = size(meta.valuemap,2);

nsee = sum(val>0,2);
valc = val(nsee==N,:);
valinc = val(nsee<N,:);

cs = [1; cumprod(card(1:end-1))];


% if nargin>2,
%     valc = [vala; valc];
% end
    
%valc = [1 1; 1 2; 2 1; 2 2; valc];
theta0 = ones(Nparam,1);
%theta0 = zeros(Nparam,1);
for i=1:Nparam,
    idx = max(meta.valuemap(:,i)-1,0)'*cs;
    vc = valc;
    vc(:,meta.valuemap(:,i)==0) = 0;
    theta0(i) = theta0(i)+sum(max(vc-1,0)*cs==idx);
end



    

cnt = 0;
cnt2 = 0;
mask = zeros(Nparam,Nparam);
for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        mask2(cnt2+1:cnt2+card(i)-1,cnt+1:cnt+card(i)) = [eye(card(i)-1) -ones(card(i)-1,1)];
        mask3(cnt2+1:cnt2+card(i)-1,cnt2+1:cnt2+card(i)-1) = 1;
        card2(cnt2+1:cnt2+card(i)-1) = card(i);
        cnt = cnt+card(i);
        cnt2 = cnt2+card(i)-1;
    end
end

Nfree = size(mask2,1);

Ntrain = size(valinc,1);




% p = ones(Nparam,1);
% k = mask*p;
% meta.p = p./k;

if nargin<5,
    meta.p = theta0./(mask*theta0);
else
    meta.p = p0;
end

go = 1;
cnt = 1;
err = 1;
while go,
    pold = meta.p;
    %hd = waitbar(0);
    theta = theta0;
%     H = zeros(Nfree,Nfree);  
    [~,m] = inferProbfast(meta,valinc);
    theta = theta+sum(m(1:Nparam,:).*m(Nnodes+(1:Nparam),:)./(ones(Nparam,1)*m(Nnodes,:)),2);
%     v = mask2*((m(1:Nparam,:)>0).*m(Nnodes+(1:Nparam),:)./(ones(Nparam,1)*m(Nnodes,:)));
%     H = H+v*v';
    if diagnostic,
        disp(sprintf('err = %f',err));
    end
%     for i=1:Ntrain,
%         waitbar(i/Ntrain,hd,sprintf('err = %f',err));
%         obs = [(1:N)' valinc(i,:)'];
%         obs = obs(obs(:,2)>0,:);
%         [~,m] = inferProb(circ,meta,obs);
%         theta = theta+m(1:Nparam).*m(Nnodes+(1:Nparam))/m(Nnodes);
%         v = mask2*((m(1:Nparam)>0).*m(Nnodes+(1:Nparam)));
%         H = H+v*v'/m(Nnodes)^2;
%     end
%     delete(hd);
    Z = mask*theta;
    meta.p = theta./Z;
    cnt = cnt+1;
    err = sqrt(mean((meta.p-pold).^2));
    if (err<1e-3)||(cnt>=50),
        go = 0;
    end
end


% for i=1:size(valc,1),
%     v = mask2*((sum((valc(i,:)'*ones(1,Nparam)).*(meta.valuemap(:,1:Nparam)>0)==meta.valuemap(:,1:Nparam),1)==N)'./meta.p);
%     H = H + v*v';
% end

id = (vala>0)*2.^(N-1:-1:0)';
id2 = (val>0)*2.^(N-1:-1:0)';
Nv = 2^N-1;

H = diag((mask2>0)*(1./meta.p))+((mask2<0)*(1./meta.p)*ones(1,Nfree)).*mask3;
H = H.*(ones(Nfree,1)*card2);
%H = diag((abs(mask2)*(1./meta.p)));
%H = zeros(Nfree,Nfree);
for i=1:Nv,
    Ni = sum(id2==i);
    if Ni>0,
        [~,m] = inferProbfast(meta,vala(id==i,:));
        v = mask2*((m(1:Nparam,:)>0).*m(Nnodes+(1:Nparam),:)./sqrt(ones(Nparam,1)*m(Nnodes,:)));
        H = H + Ni*v*v';
    end
end

pout = meta.p;
%cov = inv(eye(Nfree)+H);
cov = inv(H);
cov = mask2'*cov*mask2;

%cov = cov*(1+size(vala,1)/size(val,1));
    



        
    