function [w_out,C,Var,mc] = approx_backprop(wcond,wmarg,card,parent,wlike)

mcond = opinions_to_means(wcond); %a matrix of projected probabilities (means) of the input opinions,
                           % one column per each conditional opinion 

l= size(wcond,2)-1; % the cardinality of X
card = card(:)';
N = length(card); % N is the number of parents of the node X 

wlike = wlike(:)';
l2 = length(wlike)-1;


if l2~=l,
    error('Node likelihood cardinality incosistent with that of the conditional probability');
end

u = wmarg((1:N)'+N*card(1:N)');
S = card'./u;
wmargmean = wmarg;
wmargmean((1:N)'+N*card(1:N)') = 0;
wmargmean = wmargmean(:,1:end-1);
wmargmean = wmargmean+1./S*ones(1,max(card));

lp = card(parent);

slike = l/wlike(l+1);
mlike = wlike(1:l)+wlike(l+1)/l;

%m_out2 = prob_backward(mcond,wmargmean,card,parent,mlike);

sc = l./wcond(:,l+1); % the column with Dirichlet strengths of wcond, s=|X|/u 

Nsf = prod(card);

Ns = Nsf/lp;

wmargmean = wmargmean([1:parent-1 parent+1:N],:);
S = S([1:parent-1 parent+1:N]);

factor = [cumprod(card(2:end),'reverse') 1]; 
pfact = factor(parent);
factor = factor([1:parent-1 parent+1:N]);

card = card([1:parent-1 parent+1:N]);

N = N-1;
v1 = zeros(1,N);


M2 = zeros(Ns,Ns);
loc = zeros(lp,Ns);

mu = zeros(Ns,1);

if Ns==1,
    M2 = 1;
    loc = (1:lp)';
    mu = 1;
else

    for i=0:Ns-1,


        loc(:,i+1) = v1*factor'+pfact*(0:lp-1)'+1;

        w1 = wmargmean((1:N)'+N*v1');

        mu(i+1) = prod(w1);

        v2 = zeros(1,N); %v2 is the value of the parent to be taken in the second combination    

        for k=0:Ns-1,


            w2 = wmargmean((1:N)'+N*v2');

            mo = w1.*w2./(1+1./S)+w1.*(v1'==v2')./(S+1); %Rearanged the equation to work for dogmatic opinions where S=inf

    %             for j=1:N,
    %                 mo(j) = (wmargmean(j,v1(j)+1)*wmargmean(j,v2(j)+1)*S(j)+ wmargmean(j,v1(j)+1)*(v1(j)==v2(j)))/(S(j)+1);
    %             end
            %M2=[]; %the matrix of all the products in (61), Ns times Ns
            %matix   %you do not want to remove the previous defined
            %elements of M2

            M2(i+1,k+1) = prod(mo);

            v2(N) = v2(N)+1; % v2 is updated Ns times for each update of v1 
            for j=N:-1:2,
                if v2(j)==card(j),
                v2(j) = 0;
                v2(j-1) = v2(j-1)+1;
            else
                break;
                end
            end

        end

        v1(N) = v1(N)+1; % v1 is updated at each round 
        for j=N:-1:2,
            if v1(j)==card(j),
                v1(j) = 0;
                v1(j-1) = v1(j-1)+1;
            else
                break;
            end
        end

    end
    
end


MatLike = mlike'*mlike/(1+1/slike)+diag(mlike/(slike+1));
M = zeros(Nsf,Nsf);
for i=1:l,
    for j=1:l,
        M = M+((mcond(i,:)'*mcond(j,:)).*(1-diag(1./(sc+1)))+(i==j)*diag(mcond(i,:)'./(sc+1)))*MatLike(i,j);        
    end
end


C = zeros(lp,lp);
mc = zeros(lp,1);

for i=1:lp,
    mc(i) = mlike*mcond(:,loc(i,:))*mu;
    for j=1:i,
        C(i,j) = sum(sum(M(loc(i,:),loc(j,:)).*M2));
        if j<i,
            C(j,i) = C(i,j);
        end
    end
end




C = C - mc*mc';
% C = max(C,0);
c2 = max(diag(C),0);
C = C-diag(diag(C))+diag(c2);





sm = sum(mc);
m_out = mc/sm;


Var = max((m_out.^2*sum(C(:))-2*m_out.*sum(C,2)+diag(C))'/sm^2,0);

m_out = m_out';


s_star = sum(m_out.^2.*(1-m_out).^2)/sum(Var.*m_out.*(1-m_out))-1;

s_out = max([s_star 1./m_out]);

l = length(m_out);


w_out = [m_out-1/s_out l/s_out];