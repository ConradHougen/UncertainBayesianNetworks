function [pout,cov,cnt,mask,mask2,mask3] = em_fisher_fast(meta,val,metafisher,p0)

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

%Ntrain = size(valinc,1);




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

vala = metafisher.vala;
idvala = (vala>0)*(2.^(N-1:-1:0)');

Nc = 2^N-1;
Jc = zeros(Nfree,Nfree,Nc);
Jp = zeros(Nfree,Nfree,Nc);
jcomputed = zeros(Nc,1);

[~,m] = inferProbfast(meta,zeros(1,N));
pp = (mask2>0)*m(Nnodes+(1:Nparam));
J0 = diag(pp./((mask2>0)*m(1:Nparam)))+((pp./((mask2<0)*m(1:Nparam)))*ones(1,Nfree)).*mask3;


for i=1:N,
    loc2 = find(meta.varmap==i,1);
    v = zeros(1,N);
    v(meta.valuemap(:,loc2)>0) = 1;
    id = v*(2.^(N-1:-1:0)');
    loc = metafisher.loc{id};
    jcomputed(id) = 1;
    Jc(loc,loc,id) = J0(loc,loc);
end

idv = (val>0)*(2.^(N-1:-1:0)');
idvu = unique(idv);
h = hist(idv,0:Nc);
h = h(2:end);


for i = idvu'
    if i>0,        
        id = metafisher.dec{i};
        for j=1:length(id),
            if ~jcomputed(id(j)),
                 [~,m] = inferProbfast(meta,vala(idvala==id(j),:));
                 loc = metafisher.loc{id(j)};
                 v = zeros(Nparam,size(m,2));
                 v(loc,:) = ((m(loc,:)>0).*m(Nnodes+loc,:)./sqrt(ones(length(loc),1)*m(Nnodes,:)));
                 v = mask2*v;
                 Jc(:,:,id(j)) = v*v';   
                 jcomputed(id(j)) = 1;
            end
        end
        Jp(:,:,i) = h(i)*sum(Jc(:,:,id),3);
    end
end

Jp = sum(Jp,3);

H = diag((mask2>0)*(1./meta.p))+((mask2<0)*(1./meta.p)*ones(1,Nfree)).*mask3;
H = H.*(ones(Nfree,1)*card2);

H = H+Jp;

pout = meta.p;
%cov = inv(eye(Nfree)+H);
cov = inv(H);
cov = mask2'*cov*mask2;

%cov = cov*(1+size(vala,1)/size(val,1));
    



        
    