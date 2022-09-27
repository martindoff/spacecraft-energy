function dydt = dynamics(t,y,ft,f,gt,g, sat)
f_x = interp1(ft,f(1,:),t); % Interpolate the data set (ft,f) at time t
f_y = interp1(ft,f(2,:),t); % Interpolate the data set (ft,f) at time t
g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t


dydt(1) = y(2)  ;
dydt(2) = y(4)  ;
dydt(3) = -sat.mu*y(1)./(y(1).^2 + y(3).^2).^1.5 + 1/sat.m*(f_x.*cos(y(5)) - f_y.*sin(y(5)) ); 
dydt(4) = -sat.mu*y(3)./(y(1).^2 + y(3).^2).^1.5 + 1/sat.m*(f_x.*sin(y(5)) + f_y.*cos(y(5)) ); 
dydt(5) = y(6) ; 
dydt(6) = g/sat.J; 
dydt = dydt(:);  