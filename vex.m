function x_v = vex(A)
% Uncross operation, or ``vectorize". Input is 3 x 3 skew-sym matrix,
% output is 3 x 1 column matrix. 

x_v = [ -A(2,3); A(1,3); -A(1,2)];