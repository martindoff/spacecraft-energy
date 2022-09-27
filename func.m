function f = func(z, u, param, nu)
%% Parameters 
m = param.m; 
mu = param.mu; 
I = param.I; 
I_w = param.I_w; 
J = param.J; 
p = param.p; 
AU = 149597870; 
A = param.S_SA; 
l_SA = param.l_SA;

%% State and input unpacking
r = z(1:3,1);
v = z(4:6,1);
omega = z(7:9,1);
omega_w = z(10:12,1); 
C_ba = [z(13:15,1) z(16:18,1) z(19:21,1)]; 
r_des = z(22:24,1);
v_des = z(25:27,1);
f = u(1:3,1);
g_w = u(4:6,1); 

%% Orbital dynamics
r_sun_earth = [-AU ; 0 ; 0];
F_s = p*A*(r+r_sun_earth)/norm(r+r_sun_earth)*nu; %solar radiation force
r_dot = v; 
r_dot2 = -mu*r/norm(r)^3 + F_s + f; % EOM update %1/m*C_ba'*

%% Desired orbital dynamics 
r_dot_des = v_des; 
r_dot2_des = -mu*r_des/norm(r_des)^3;

%% Attitude dynamics
r_b = C_ba * r; 
ratio_p = 161/85670; % ratio time in half penumbra vs orbit periode 
arc_p = 2*pi*norm(r)*ratio_p; % arc length of half penumbra section 
dec = l_SA / arc_p;
if (nu < 1) % only in penumbra 
tau_SA = [0 ; 0; sign(r(2))*norm(F_s)*l_SA^2/6*dec]; % solar torque
else 
tau_SA = [0 ; 0; 0];
end 
tau_grav= (3*mu/norm(r)^5)* cross(r_b) * J * r_b ; % gravity gradient
tau_ext = tau_SA + tau_grav; % external torques 
omega_dot = I^(-1)*( - cross(omega) *( J * omega + I_w*omega_w )   - g_w + tau_ext); % EOM update


%% Wheel dynamics
 omega_w_dot  =  I_w^(-1)*g_w - omega_dot ; % EOM update

%% Rotational kinematics
C_ba_dot = - cross(omega)*C_ba; % EOM update

%% State repacking 
f = [r_dot; r_dot2; omega_dot; omega_w_dot; C_ba_dot(:,1); C_ba_dot(:,2); C_ba_dot(:,3); r_dot_des; r_dot2_des]; 

end

