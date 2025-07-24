function node = cdh_theta_mle(node, val, meta)
    
N = length(node);

count = zeros(1, length(meta.varmap));

for i=1:N
    valmapidx = find(meta.varmap == i);
    np = length(node(i).parents);
    
    % If no parents, just look at values that node i takes on
    if np == 0
        for j=valmapidx   
            count(j) = sum(val(:,i) == meta.valuemap(i,j));
        end
        node(i).mle(1:node(i).cardinality) = count(valmapidx)/sum(count(valmapidx));
    % If parents
    else
        card = node(i).cardinality;
        for parent=node(i).parents
            for k=1:length(valmapidx)/card
                for j=valmapidx((k-1)*card+1:k*card)
                    % validx is the indices where the parent takes on a
                    % given value, in the val array. Since we are
                    % conditioned on the parent being equal to some value,
                    % we only need to look at the subset of rows of val 
                    % which match this condition
                    validx = val(:,parent) == meta.valuemap(parent,j);
                    % Looking only at the conditioned rows, count the
                    % number of times that node i takes on value in
                    % valuemap(i,j)
                    count(j) = sum(val(validx,i) == meta.valuemap(i,j));
                    node(i).mle(j == valmapidx) = count(j)/length(val(validx,i));
                end         
            end
        end
                    
                
    end    
    
end

end