% Function to compute sigmoid(x)
function [sigm] = cdh_sigmoid(x)

sigm = exp(x)./(1+exp(x));

end