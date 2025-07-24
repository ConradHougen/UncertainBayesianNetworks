% File for simply generating and saving GTBNs and val matrices.
function [out] = cdh_generate_gtbns_vals(rg,file,Nruns,Ntrain,Fobs)
    
    node = makenetwork(file);
    Nvars = length(node);
    
    out.Nvars = Nvars;
    out.val = zeros(Ntrain, Nvars, Nruns);
    out.obsmask = zeros(Ntrain,Nvars,Nruns);
    
    for i = 1:Nruns
        node = createBN(node,rg);
        out.node(i,:) = (node);
        
        [out.circ, meta] = cdh_net2circuit(node);
        out.meta(i) = setup(out.circ,meta);
        
        out.val(:,:,i) = cdh_sim_observations(node, Ntrain);
        out.obsmask(:,:,i) = rand(size(out.val(:,:,i))) < Fobs;
        %out.obsmask(1:20,:,i) = repmat([1 1 1 1 1 1 1 1 1], 20, 1);
        %out.obsmask(21:end,:,i) = repmat([1 0 0 1 0 0 1 0 0], Ntrain-20, 1);
        out.val(:,:,i) = out.val(:,:,i) .* out.obsmask(:,:,i);
    end
    
    [out.fpidx,out.ifpidx] = cdh_find_free_params(node, out.meta(1));
    
end

