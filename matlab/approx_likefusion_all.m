function [wout,wbeta] = approx_likefusion_all(win)

% win is a matrix of opinions on the same variable X, one row per each
% opinion

l=length(win(1,:))-1; %l is the cardinality of X
c=length(win(:, 1)); %the number of opinions to be fused

%epsilon = 2; 

%a=ones(l,1).*1/l; %the base rate of X

%wintransposed=win';

%bin=wintransposed(1:l, :); %the belief part of win columnwise, one column for each (belief) opinion

%uin=win(:,l+1); %the column with the uncertainties of win

%uinarep=repmat(uin.*a, 1, c); %replicate uin.*a to a matrix with c same columns

%min = bin+uinarep; %a matrix of projected probabilities (means) of the input opinions,
                   % one column for each opinion i
                   
[min,uin] = opinions_to_means(win);                   
                   
mprod = prod(min, 2);  %the product of the mean values at a single element
mout = mprod/sum(mprod); %the mean distribution of the output opinion 

%sin = (1/l)./uin; %the column of initial Dirichlet strengths s^i_X
sin = l./uin;  %LMK - Note that S = K/u
moutrep = repmat(mout, 1, c);

sumelement = sum((moutrep.^2)./min,1);%calculating the sum in the parentheses in step 3 operator 1 equation
                                  % sumelement is a column of 
                                  %LMK added the second argument in the sum
                                  %to ensure always a sum over columns
% minfirstrow = min(1,:); %we choose one (the first) value of X to work with
% 
% moutfirstelement = mout(1,:);

%sdesired = ((1/(((1-2*moutfirstelement)./minfirstrow)+sumelement)*sin))*moutfirstelement*(1-moutfirstelement))-1;
%Step 3 in Operator 1 - to be checked!

%sdesired = 1/sum(((1-2*moutfirstelement)./minfirstrow+sumelement)*moutfirstelement/(1-moutfirstelement)./(sin+1),2)-1,

%LMK - determine sdesired to fit the variance for each of the c element values of the variable
mout2 = mout*ones(1,c);
sdesired = 1./sum(((1-2*mout2)./min+ones(l,1)*sumelement).*mout2./(1-mout2).*(1./(ones(l,1)*sin+1)),2)-1;

if nargout>=2,
    sdesired = max([sdesired 1./mout 1./(1-mout)],[],2);
    mat = sdesired*mout';
    alpha = [diag(mat) diag(mat*(1-eye(l)))];
    wbeta = [alpha-1 2*ones(l,1)]./(sum(alpha,2)*[1 1 1]);
end



if isinf(max(sdesired)),
    sactual = inf;
else
    
    %LMK - Least squares fit of the approximated variance and the fitted Dirichlet variance
    sdesired = sum(mout.^2.*(1-mout).^2,1)/sum(mout.^2.*(1-mout).^2./(sdesired+1),1)-1;

    sactual = max(max(sdesired, 1./mout)); %step 4, operator 1
end

%alphaout = (mout.*sactual)'; %step 5, operator 1


wout = [mout'-1/sactual l/sactual];

     

