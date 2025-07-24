function we_so = cdh_gen_subjective_opinions(mtheta, rtheta)

    [M,N] = size(mtheta);
    we_so = zeros(M*N, 3);
    
    % For each observation
    for j = 1:N
        % First, compute Dirichlet strength
        m = mtheta(:,j);
        v = rtheta(:,j);
        s = (m.*(1-m))./v - 1;
        
        % Compute subjective opinions
        we_so((j-1)*M+1:j*M, :) = [m-1./s, 1-m-1./s, 2./s];
    end
    

end