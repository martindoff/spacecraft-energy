function u   = control(z, param)

%% Unpacking 
r = z(1:3,1);
r_dot = z(4:6,1);
omega = z(7:9,1);
C_ba = [z(13:15,1) z(16:18,1) z(19:21,1)]; 
r_des = z(22:24,1);
r_des_dot = z(25:27,1);
m=param.m;
AU = param.AU;
 

%% Orbital control

% Gains
kp=.1; 
kd=.1; 

% Control
f = m*C_ba*(kp*(r_des - r) + kd*(r_des_dot - r_dot)) ;

%% Attitude control

% Gains 
k = 0.012/10;
u_max = 3*1200e-3;
w_max = u_max*(1-k);
Kp = k*u_max/pi;

% Desired 
r_sun_earth = [-AU ; 0 ; 0];
x_a = [1;0;0]; y_a =[0;1;0]; % inertial basis vectors
s = -(r + r_sun_earth)/norm(r + r_sun_earth); % unit vector pointing towards sun
s_y = [-s(2) ; s(1); 0];

C_da = [s'*x_a     s'*y_a      0;...
        s_y'*x_a   s_y'*y_a    0;...
             0          0      1 ];
         
% Control 
C_bd  = C_ba*C_da'; % error DCM
P = .5*(C_bd - C_bd'); % projection from SO(3) -> so(3) 
g = -Kp*vex(P)/sqrt(1+trace(C_bd))-sigma(omega, w_max); % control torque


u = [f;g]; 
end



