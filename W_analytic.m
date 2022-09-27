function W = W_analytic(x)
%W_SIMPLE Summary of this function goes here
%   Detailed explanation goes here
eps = 0.4586887; 
W = (1+eps)*log(6/5*x./log( 12/5*x./log(1+12/5*x) )) - eps*log(2*x./log(1+2*x)); 



end

