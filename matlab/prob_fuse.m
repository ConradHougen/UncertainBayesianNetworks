function m_out = prob_fuse(m_in)

% the operation of fusion for probability distributions, 
% m_in - matrix with one distribution on the same variable per each row 

a = prod(m_in,1); % a row with products of the columns of m_in, 
                % one product for each value of the variable
%b = prod(1-m_in);
m_out = a/sum(a); % a normalization of the row a, gives a single probability distribution 
                  % over the same variable as the initial ones 