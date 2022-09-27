%% Initialise all parameters 

%% Physical parameters
mu=3.986004418e+14*1e-9;  
q_e = 1.602176634e-19 ; % electric charge of an electron (C)
sigma = 5.670e-8;
G_0 = 63300000.0; 
R_sun = 696000000.0;
R_earth = 6378000.0; 
AU = 150000000000.0; 
vehicle.AU = AU; 

%% Vehicle dynamics 
I=diag([1.1756e+04; 1.1756e+04; 5100]);
I_w=diag([22.5000; 22.5000; 22.5000]); 
J = I + I_w; 
m=3551.0;
vehicle.I = I; 
vehicle.I_w = I_w;
vehicle.J = J; 
vehicle.m = m; 
vehicle.mu = mu;
omega0 = 0; 
z=42164000.0; 
r0 = [-42164 ; 0 ; 0];
v0 = [0; 3.0690; 0];
omega0 = [0;0;0];
omega_w0 = [0; 0; 0]; 
C_ba0 = [1; 0; 0; 0; 1; 0; 0; 0; 1]; 

%% Battery 
R_bat= 120e-3*28/5/10;  
U_bat = 28*3.6; 
C_th_bat= 174581.8;
T_max_bat= 15.0 + 273.15;
T_min_bat = 5 +273.15;
T_0_bat = 11.73 + 273.15;
Q0=3600*250;
G1_th = 913.42/2/10*5;
G2_th = 205*0.05/0.02*5; % heat pipe (aluminum)  
T_fluid = 273.15 +0; 
G_th = G1_th + G2_th; 
E_max = 830822400/150*10*2*2/1.5; 
E_min = 0.2*E_max; 
eff=0.5;
P_comp = 3500; % note : in paper, bat -> comp 


%% External surface 
C_th_surface = 2760453.0/50; 
epsilon_surface = 0.81; 
S_surface = 17.1; 
phi=0; 
gamma_surface = 0.87;  
th_surface.d = d; 
th_surface.C_th = C_th_surface; 
th_surface.epsilon = epsilon_surface;
th_surface.sigma = sigma; 
th_surface.S = S_surface; 
th_surface.G_0 = G_0; 
th_surface.R_sun = R_sun;
th_surface.AU = AU; 
th_surface.phi = phi; 
th_surface.R_earth = R_earth; 
th_surface.z=z; 
th_surface.gamma = gamma_surface;
th_surface.G1_th = G1_th; 
th_surface.T_mean_bat = (T_max_bat + T_min_bat)/2; 

%% Solar array 
T_0_SA = -39.04+273.15; 
gamma_SA = 0.5; 
C_th_SA= 87000.0; 
epsilon_SA=0.85;
S_SA= 108.0; 
l_SA = 2.0*18.0 + 3; %length lever arm for solar force
p = 4.5*10^(-6); %solar pressure
vehicle.p = p; 
vehicle.S_SA = S_SA;
vehicle.l_SA = l_SA; 
th_SA.d = d; 
th_SA.C_th = C_th_SA; 
th_SA.epsilon = 2*epsilon_SA;
th_SA.sigma = sigma; 
th_SA.S = S_SA; 
th_SA.G_0 = G_0; 
th_SA.R_sun = R_sun;
th_SA.AU = AU; 
th_SA.phi = phi; 
th_SA.R_earth = R_earth; 
th_SA.z=z; 
th_SA.gamma = gamma_SA;
param.V_oc_cell_SA =2.665;
param.T_ground_SA =  301.15;
param.alpha_V_oc_cell_SA = -0.0062;
param.I_sc_cell_SA = 0.5456; 
param.T_ground_SA = 301.15;
param.alpha_I_sc_cell_SA = 0.00032; 
param.k_B = 1.380649e-23; 
param.q_e =  1.602176634e-19; 
param.N_p_cell_SA = 450; 
param.N_s_cell_SA = 75; 
param.I_0 = 0.0006574584908885646/1.5; 
param.G_0 = G_0; 
param.R_sun = R_sun; 
 
%% Propulsion and actuation
Kv = 3600; % motor constant
k_a = 60/(2*pi*Kv)/10; % motor torque constant 
R_a = 1; % armature resistance 
k_w = k_a; % back emf constant
V_HET = 100.0; % operating voltage HET (V)
eta_HET = 0.6 ; % efficiency of HET (-)
N_HET = 6.0 ; % number of thrusters (-)
q_HET = (0.8*1 + 0.2*2)*q_e ; %average charge of Xenon atoms (C) 
M_HET = 131.293*1.660538921e-27 ; %atomic mass Xenon (kg)


%% Solver 
gamma=10*0; 
P=1e+4; P_max=U_bat^2/(4*R_bat); E_=E_max; % I=100; I_max=U_bat/(2*R_bat); T_=292;