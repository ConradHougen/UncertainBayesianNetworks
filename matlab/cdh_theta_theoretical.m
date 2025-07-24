% Compute the theoretical mean/covariance for the theta variables, based on
% Dirichlet posterior
function [thetaBar, thetaCov] = cdh_theta_theoretical(node, val, meta)
    
Nvar = length(node);
Ntheta = length(meta.varmap);
count = zeros(1, Ntheta);
alpha = zeros(1, Ntheta);

thetaBar = zeros(Ntheta, 1);
thetaCov = zeros(Ntheta);

% Remove incomplete rows of val
val = val(~any(val == 0, 2),:);

for i=1:Nvar
    valmapidx = find(meta.varmap == i);
    np = length(node(i).parents);
    
    % If no parents, just look at values that node i takes on
    if np == 0
        for j=valmapidx   
            count(j) = sum(val(:,i) == meta.valuemap(i,j));
        end
        alpha(valmapidx) = count(valmapidx)+1;
        a0 = sum(alpha(valmapidx));
        atilde = alpha(valmapidx)/a0;
        
        % Mean of Dirichlet
        thetaBar(valmapidx) = alpha(valmapidx)/a0;
        
        % Covariance of Dirichlet
        for m = 1:length(valmapidx)
            for n = 1:length(valmapidx)
                if m == n
                    thetaCov(valmapidx(m), valmapidx(n)) = atilde(m)*(1-atilde(m))/(a0 + 1);
                else
                    thetaCov(valmapidx(m), valmapidx(n)) = -alpha(valmapidx(m))*alpha(valmapidx(n))/(a0^2 *(a0+1));
                end
            end
        end
        
    % If parents
    else
        card = node(i).cardinality;
        for parent=node(i).parents
            for k=1:length(valmapidx)/card
                idxarr = valmapidx((k-1)*card+1:k*card);
                for j=idxarr
                    % validx are the indices where the parent takes on a
                    % given value, in the val array. Since we are
                    % conditioned on the parent being equal to some value,
                    % we only need to look at the subset of rows of val 
                    % which match this condition
                    validx = val(:,parent) == meta.valuemap(parent,j);
                    
                    % Looking only at the conditioned rows, count the
                    % number of times that node i takes on value in
                    % valuemap(i,j)
                    count(j) = sum(val(validx,i) == meta.valuemap(i,j));                   
                end
                
                alpha(idxarr) = count(idxarr) + 1;
                a0 = sum(alpha(idxarr));
                atilde = alpha(idxarr)/a0;
                
                % Mean of Dirichlet
                thetaBar(idxarr) = alpha(idxarr)/a0;
                
                % Covariance of Dirichlet
                for m = 1:length(idxarr)
                    for n = 1:length(idxarr)
                        if m == n
                            thetaCov(idxarr(m), idxarr(n)) = atilde(m)*(1-atilde(m))/(a0 + 1);
                        else
                            thetaCov(idxarr(m), idxarr(n)) = -alpha(idxarr(m))*alpha(idxarr(n))/(a0^2 *(a0+1));
                        end
                    end
                end
            end
        end
                    
                
    end    
    
end

end