function [p,cov,hits,alpha] = ep(meta,val)

maxcnt = 20;
cnt2 = 1;

Ntrain = size(val,1);
N = size(val,2);
nparam = length(meta.varmap);
nnodes = size(meta.valuemap,2);

v = reshape(meta.valuemap(:,1:nparam),[N 1 nparam]); 
v = repmat(v,[1 Ntrain 1]);

na = sum(meta.valuemap(:,1:nparam)>0,1);
na = reshape(na,[1 1 nparam]);
na = repmat(na,[1 Ntrain 1]);

valn = val;
valn(val==0) = NaN;
hits = sum(repmat(valn',[1 1 nparam])==v,1)==na;
hits = double(reshape(hits,[Ntrain nparam]));

alpha = sum(hits,1)'+1;

cnt = 0;
mask = zeros(nparam,nparam);

card = max(meta.valuemap(:,1:nparam),[],2);


for i=1:N,
    loc = find(meta.varmap==i);
    Ne = length(loc)/card(i);
    for j=1:Ne,
        mask(cnt+1:cnt+card(i),cnt+1:cnt+card(i)) = 1;
        cnt = cnt+card(i);
    end
end

M = ones(nparam,1)*meta.varmap;
maskfull = M==M';


M1 = eye(nparam);
M2 = mask-M1;
M3 = maskfull-mask;


go = 1;


while go,

    alpha_old = alpha;

    for i=1:Ntrain,

        alpham = alpha-hits(i,:)';
        S = mask*alpham;
        p = alpham./S;

        meta.p = p;
        [~,m] = inferProbfast(meta,val(i,:));
        m2 = m(1:nparam).*m((1:nparam)+nnodes);

        p1 = p*ones(1,nparam);
        p2 = ((alpham+1)./(S+1))*ones(1,nparam);
        p22 = (alpham./(S+1))*ones(1,nparam);
        p3 = ((alpham+2)./(S+2))*ones(1,nparam);
        p32 = ((alpham+1)./(S+2))*ones(1,nparam);

        mn = ((p2.*M1+p22.*M2+p1.*M3)*m2)/m(nnodes);
        vv = ((p3.*p2.*M1+p32.*p22.*M2+p2.*p1.*M3)*m2)/m(nnodes);

        S = (mask*(mn.*(1-mn).*(mn-vv)))./(mask*(mn.*(1-mn).*(vv-mn.^2)));
        alpha = mn.*S;
        
       
       % adjust = mask*max(max(hits(setdiff(1:Ntrain,i),:))'-alpha+1,zeros(nparam,1));
        adjust = 0;
        alpha = alpha+adjust; %prevent negative alpahm
        
        if sum(adjust)>0, %sum(sum(alpha*ones(1,Ntrain)-hits'<0,1),2)>0,
            %disp('adjust');
        end
        
        hits(i,:) = (alpha-alpham)';
        

        p = mn;

    end
    
    %[alpha'; alpha_old'],
    
    
    delta = norm(alpha-alpha_old);
    go = (delta>1e-4)&&(cnt2<maxcnt);
    cnt2 = cnt2+1;

    
end

disp(num2str(cnt2));


cov = (diag(p./(S+1))-(p./(S+1))*p').*mask;



