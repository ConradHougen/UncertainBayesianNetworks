function mout = prob_backward(mc,m_in,card,parent,ml)

% the function determines the lambda-message that X_i sends to its parent node X_k 

% m_in - the matrix with the \pi-messages of the (other) parents of X,
% one row per each setting of the parents

% mc = node(i).p - the initial matrix of conditional prob. of X_i given the
% parents, one column per parent

% card is the cardinality of each parent variable

% parent - the index k of the parent the node X_i is sending the message

% ml- the matrix of the \lambda values of X_i

% prob_backward(node(id).p,node(id).parentmsg,sendid,node(id).back)

n = size(mc,1);

if nargin<4,
    ml = 1/n*ones(1,n);
end

N = length(card);
Ns = prod(card);

loop = setdiff(1:N,parent);
v = zeros(1,N);
loc = zeros(Ns,1);
mo = ones(Ns,1);
for i=0:Ns-1,
    loc(i+1) = v(parent)+1;
    for j=loop,
        mo(i+1) = mo(i+1)*m_in(j,v(j)+1);
    end
    v(N) = v(N)+1;
    for j=N:-1:2,
        if v(j)==card(j),
            v(j) = 0;
            v(j-1) = v(j-1)+1;
        else
            break;
        end
    end
end

mc = ml*mc;

N = card(parent);
mzx = zeros(1,N);
for i=1:N,
    mzx(i) = mc(loc==i)*mo(loc==i);
end
mout = mzx; %/sum(mzx);

% N = length(m_in);
% Ns = 2^(N+1);
% 
% mo = ones(1,Ns);
% for i=[1:parent-1 parent+1:N],
%     ix = bitand(bitshift(0:Ns-1,i-N),1);
%     mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
% end
% ix = bitand(bitshift(0:Ns-1,-N),1);
% iy = bitand(0:Ns-1,Ns/2-1)+1;
% m = [1-ml ml];
% mc = mc';
% mo = mo.*((1-mc(iy)).^(1-ix)).*(mc(iy).^ix).*m(ix+1);
% 
% ix = bitand(bitshift(0:Ns-1,parent-N),1);
% 
% mzxo = sum(mo(ix==1));
% mznxo = sum(mo(ix==0));
% 
% 
% mout = mzxo/(mzxo+mznxo);


