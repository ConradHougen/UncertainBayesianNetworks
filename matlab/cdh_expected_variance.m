function evavg = cdh_expected_variance(w)

    v = [1; 0; 0.5];
    m = w * v;
    s = 2 ./ w(:,3);
    
    ev = (m.*(1-m))./(s+1);
    evavg = mean(ev);
    
end