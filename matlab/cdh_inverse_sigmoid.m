function sigminv = cdh_inverse_sigmoid(y)

sigminv = log(y./(1-y));

end