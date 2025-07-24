function [ m_in, u_in ] = opinions_to_means( w_in )

%for a matrix w_in of opinions on X given columnwise, this function 
%returns the matrix of the corresponding mean distributions m_in

l=length(w_in(1,:))-1; %l is the cardinality of X
c=length(w_in(:,1)); %the number of conditional opinions in the input


a=ones(l,1).*1/l; %the base rate of X

w_intransposed=w_in';

b_in=w_intransposed(1:l, :); %the belief part of w_in columnwise, 
                             %one column per each (belief) opinion

u_in=w_in(:,l+1); %the column with the uncertainties of w_in

%u_in_rep=repmat(u_in.*a, 1, c); %replicate u_in.*a to a matrix with c same columns
u_in_rep = a*u_in';  %LMK - rhe uncertainty is a function of opinion (column) and the prior is a function of class (row);

m_in = b_in + u_in_rep; %a matrix of projected probabilities (means) of the input opinions,
                           % one column per each conditional opinion 
                           
u_in = u_in';  %LMK- to be consistent with the columns representing opinions


end

