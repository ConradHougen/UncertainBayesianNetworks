function [fpidx, ifpidx] = cdh_find_free_params(node, meta)

    Nparams = length(meta.varmap);
    
    % Free parameters are all parameters except for the final parameter in
    % each set, so find those locations and invert them
    inversefpidx = false(1,Nparams);
    
    i = 1;
    while i <= Nparams
        % Determine which BN variable param i belongs to
        var = meta.varmap(i);
        % Determine cardinality of that BN variable
        card = node(var).cardinality;
        % Add to inverse of free params
        i = i + card - 1;
        inversefpidx(i) = true;
        
        % Go to the next parameter
        i = i + 1;
    end
    
    fpidx = find(xor(true(1,Nparams), inversefpidx) == true);
    ifpidx = [0 find(inversefpidx == true)];
    
end