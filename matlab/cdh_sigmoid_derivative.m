function dsigm = cdh_sigmoid_derivative(x)

dsigm = exp(x)./(1 + exp(x)).^2;

end