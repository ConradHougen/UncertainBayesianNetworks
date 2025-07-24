function wout = fixw(w)

m = w*[eye(2); 0.5 0.5];
s = 2./w(:,3);
s = max([1./m s],[],2);
wout = [m-1./s*[1 1] 2./s];