% Apply probability constraints to theta (KKT conditions)
function theta = cdh_apply_constraints(theta, meta, node)

N = length(node);

for k = 1:N
    % Get indices in theta associated with node k
    loc = meta.varmap == k;
    % Get cardinality of node k and reshape the associated part of theta
    card = node(k).cardinality;
    if card <= 0
        fprintf(1, 'Error: cardinality %d invalid', card);
        return;
    end
    
    th = reshape(theta(loc), card, []);

    % Each column of th should now be a probability distribution
    for j = 1:size(th,2)
        % First, add "lambda" to make all elements nonnegative
        [~,lambidx] = min(th(:,j));
        lambidx = lambidx(1); % In case of multiple minimum elements
        if th(lambidx, j) < 0
            th(:,j) = th(:,j) - th(lambidx, j);
        end
        
        % Now, normalize the rest of the elements to a probability
        % distribution
        s = sum(th(:,j));
        % Handle corner case where entire column is zero
        if s < 1e-6
            th(:,j) = ones(size(th(:,j))) / card;
        end
        th(:,j) = th(:,j) / s;
        
    end

    theta(loc) = reshape(th, [length(theta(loc)), 1]);
    
end
    
end