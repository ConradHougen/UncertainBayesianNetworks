function th = cdh_kkt_constrain(th, ifpidx)

    for i = 1:length(ifpidx)-1
        fpset = ifpidx(i)+1:ifpidx(i+1)-1; % Indices of free params in set i

        thfp = th(fpset);
        minth = min(thfp, 0);

        nfp = length(fpset);
        for j = 1:nfp
            if min(minth) < 0
                % Zero out the min params
                thfp = thfp - minth;
            end

            % Confirm that other ineq constraint is satisfied
            if sum(thfp) <= 1
                break;
            else
                % If sum of free params is greater than 1, need to reduce
                % the inactive constraint variables by that magnitude
                loc = find(minth == 0);
                mu = (sum(thfp) - 1)/length(loc);
                thfp(loc) = thfp(loc) - mu;
            end
            
            minth = min(thfp, 0);
        end
        
        % Set non-free param as 1 minus free params
        th(fpset) = thfp;
        th(ifpidx(i+1)) = 1 - sum(thfp);
    end

end
