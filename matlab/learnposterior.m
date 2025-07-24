function [m, r] = learnposterior(val)

    % This function only runs for a 2-variable Bayes Net!
    [Ntrain,Nvars] = size(val);
    if Nvars ~= 2
        fprintf(1, 'Warning: %d variables, so learnposterior will not run\n', Nvars);
        m = 0;
        r = 0;
        return;
    end

    % Generate nin from val
    zeroX = find(val(:,1) == 0); % zero means unobserved
    zeroY = find(val(:,2) == 0);
    zeroXorY = union(zeroX, zeroY); % All rows with at least one zero
    Xobs = setdiff(1:Ntrain, zeroX);
    Yobs = setdiff(1:Ntrain, zeroY);
    XYobs = setdiff(1:Ntrain, zeroXorY);
    Xtrue = find(val(:,1) == 1);
    Xfalse = find(val(:,1) == 2);
    Ytrue = find(val(:,2) == 1);
    
    nin.Na = length(XYobs);
    nin.Tx = length(Xobs);
    nin.Ny = length(intersect(zeroX, Yobs));
    nin.Tyx = length(intersect(Xtrue, Yobs));
    nin.Tynx = length(intersect(Xfalse, Yobs));
    nin.nx = length(Xtrue);
    nin.nyx = length(intersect(Xtrue, Ytrue));
    nin.nynx = length(intersect(Xfalse, Ytrue));
    nin.ny = length(intersect(zeroX, Ytrue));
    
    [m, r, aout] = learnpostv2(nin);
    
end