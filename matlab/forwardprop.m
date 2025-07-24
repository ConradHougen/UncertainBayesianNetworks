function [w_out,M,M2] = forwardprop(wcond,wmarg,card)

%LMK updated 11/20/19 at 3:16pm

%wcond is a matrix of conditional opinions, one per each row (one per each setting of the
%parents of X)

%wmarg is a matrix of the marginal opinions of the parents of X given the evidence
%a row of this matrix looks like:b_{u_i=0|o_i},...,b_{u_i=k_i|o_i},u_{u_i|o_i}
%the 'empty' spaces are zeros

%card is the cardinality of each parent variable, a row of numbers?

%wout is a column with the opinion of the variable X given the evidence


mcond = opinions_to_means(wcond); %a matrix of projected probabilities (means) of the input opinions,
                           % one column per each conditional opinion 

l= size(wcond,2)-1; % the cardinality of X
N = length(card); % N is the number of parents of the node X 

%The allocation of u and S as empty is not needed if we remove the loop

%u=[]; % a matrix of uncertanties of each parent 
%S=[]; % a matrix of Dirichlet strengths for each parent 

%If not replacing the loop below, the wmargmean should be preallocated
%wmargmean=zeros(N, max(card)); %matrix of the conditional opinions' means

%replacing the loop below 
u = wmarg((1:N)'+N*card(1:N)');
S = card'./u;
wmargmean = wmarg;
wmargmean((1:N)'+N*card(1:N)') = 0;
wmargmean = wmargmean(:,1:end-1);
wmargmean = wmargmean+1./S*ones(1,max(card)); %.*(wmargmean>0);  Got rid of the extra multiplication as the values that are beyond the cardinality are not used! %Not sure if the multiplication with (wmargmean>0) is necessary as the extra elements should be ignored


% for i=1:N,
%     u(i) = wmarg(i, card(i)+1);
%     S(i) = card(i)/u(i); %Dirichlet strength of the i-th parent opinion 
%     for j=1:card(i),
%         wmargmean(i,j)=wmargmean(i,j)+wmarg(i,j)+wmarg(i, card(i)+1)/card(i);  
%     end            
% end


m_out = prob_forward(mcond,wmargmean,card); %step 2, operator 2


sc = l./wcond(:,l+1); % the column with Dirichlet strengths of wcond, s=|X|/u 


%the following code creates a matrix M2, dimension Ns times Ns, 
% where each M2(i,k) is one product in (61):

if N==0, % if X is a root node
    w_out = wcond';
else
    Ns = prod(card); % Ns is the number of conditional opinions at X

    v1 = zeros(1,N); %v1 is the value of the parent to be taken in the first combination 
               
    M2 = zeros(Ns,Ns); %[]; %Better to define the full size of the array instead of having to allocate memory inside the loops
    
    for i=0:Ns-1,
        
        v2 = zeros(1,N); %v2 is the value of the parent to be taken in the second combination  
        
        w1 = wmargmean((1:N)'+N*v1');
        
        for k=0:Ns-1,
             
            %mo= zeros(1,N); %[]; %the column of (64) for each parent  %Better to define the full size of the array instead of letting it grow in the loop
            
            % It hink this helps to remove the loop below to define mo for each
            % element
            
            
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
    
       
%the following code takes a row in the mcond matrix (a particular x) and
%calculates the corresponding variance:

    %mcondactive=[]; %No need to preallocate as you are defining the
    %complete variable below

    Var = zeros(1,l); %m_out is a row vector and thus make Var a row vector
    for i=1:l,
        %mcond is akin to the tranpose of pcond
        mcondactive = mcond(i, :)'; %takes the ith row corresponding to ith value of X
        % to calculate the variance and Dirichlet strength
        % of the output
        M = mcondactive*mcondactive'+diag(mcondactive.*(1-mcondactive)./(sc+1)); % equation 62 
        %M has dimension Ns times Ns

        EM = M.*M2;
        EM = sum(EM(:));%equation (61)

        Var(i) = EM - m_out(i)^2; %the variance, Step 4, Operator 2
    end

    s_star = sum(m_out.^2.*(1-m_out).^2)/sum(Var.*m_out.*(1-m_out))-1;
    
        %s_star = m_out(i)(1-m_out(i))/Var - 1;

    s_out = max([s_star 1./m_out]);

   % end


   %alpha_out = m_out *s_out; %Step 6  %Note s_out is a scalar and thus no need for a pointwise multiplication
   %The appove line is really not needed to define w_out below.

    w_out = [m_out-1/s_out l/s_out];
    
end





    
    