function rmse = cdh_rms_error(p, w)

    v = [1; 0; 0.5];
    
    ppred = w * v;
    m = length(ppred);
    
    rmse = sqrt(sum((p-ppred).^2)/m);   

end