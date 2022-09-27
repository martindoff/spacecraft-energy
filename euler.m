function eul = euler(C)
% DCM to Euler
phi = atan2(C(1,2),C(1,1));
theta = asin(-C(1,3)); 
psi = atan2(C(2,3),C(3,3)); 
eul = [phi; theta; psi]; 
end 
