%% Optimisation - CVX (Single-shot)
E_max = 830822400/150*10*2*2;
E_min = 0.3*E_max; 
E_0 = E_min;
T_0 = T_0_bat;

%% Resize problem
N = L; 
Phi=ones(N,1);
a_bat=(1 - exp(-G_th*t/C_th_bat))/G_th ;
b_bat = T_0*exp(-G_th*t/C_th_bat);


%% CVX solver 
cvx_begin quiet
    %cvx_precision high
    cvx_solver mosek
    variables P_out(N) P_load(N) E(N) P_heat(N) P_m(N) T_bat(N) P_b(N) xi 
    minimize( sum(-P_load*P )*d + gamma*(E_min -xi))
    subject to
           E*E_ <=   E_max*Phi 
           E_min * Phi <= E*E_
           E(1)*E_ == E_0
           %E(1)*E_ == E(end)*E_
           xi==E(end)*E_
           P_out*P + V_d.*I_s  ==  P_load*P + P_a + P_HET  +  P_heat*P + P_comp/2  %I_L - I_ds*(exp(V_x)-1.0) + abs(P_heat*P)/U_heat
           P_load*P <= P_max
           P_load*P >= 0.0
           P_heat*P <= P_max
           P_heat*P >= 0
           P_out*P <= P_max
           %P_out*P >= -P_max
           %P_b*P <= U_bat^2/(2*R_bat)
           %P_b*P >= -U_bat^2/(2*R_bat)
           P_m*P >= max(P_b*P, P_b*P*eff)
           P_b*P >= U_bat^2/(2*R_bat)*(1 -sqrt(1- 4*R_bat/U_bat^2*P_out*P)) 
           E(2:end)*E_ == E(1:end-1)*E_ -d*P_m(1:end-1)*P
           %P_out(end)*P == 0
           %P_b(end)*P == 0
           a_bat.*(P_comp+P_heat*P + G1_th*T_surface + G2_th*T_fluid) + b_bat <= Phi*(T_max_bat)
           T_min_bat*Phi <= a_bat.*( P_comp+ P_heat*P + G1_th*T_surface + G2_th*T_fluid) + b_bat
 
cvx_end
cvx_status

% Rescaling
T_bat = a_bat.*(P_comp+P_heat*P + G1_th*T_surface + G2_th*T_fluid) + b_bat;
P_out = P_out*P; P_load = P_load*P;  P_heat = P_heat*P;   %I_heat = I_heat*I; I_out = I_out*I;
P_m = P_m*P; E= E_*E; P_b = P_b*P;
P_s = V_d.*I_s;


%% Plots: 
figure()
subplot(3,2,1)
plot(t/3600, P_s/1000, '-b', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$P_s$ (kW)', 'fontsize',15,'Interpreter','latex')
%title('Solar array power', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 
subplot(3,2,2)
plot(t/3600, T_bat-273.15, '-b', 'LineWidth', 2)
hold on 
plot(t/3600, ones(size(t))*(T_max_bat-273.15), '--r', 'LineWidth', 2)
plot(t/3600, ones(size(t))*(T_min_bat-273.15), '--r', 'LineWidth', 2)
axis([0 t(end)/3600 (T_min_bat -273.15)*0.8 (T_max_bat-273.15)*1.2])
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$T_c$ ($^\circ$C)', 'fontsize',15,'Interpreter','latex')
%title('Electronics temperature', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 
subplot(3,2,3)
plot(t/3600, E/1e+6, '-b', 'LineWidth', 2)
hold on 
plot(t/3600, ones(size(t))*E_max/1e+6, '--r', 'LineWidth', 2)
plot(t/3600, ones(size(t))*E_min/1e+6, '--r', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$E$ (MJ)', 'fontsize',15,'Interpreter','latex')
%title('Battery energy', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
axis([0 t(end)/3600 E_min/1e+6*0.8 E_max/1e+6*1.2])
grid on 
subplot(3,2,5)
plot(t(1:end-1)/3600, P_out(1:end-1)/1000, '-b', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$P_{out}$ (kW)', 'fontsize',15,'Interpreter','latex')
%title('Battery output power', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 
subplot(3,2,4)
plot(t(1:end-1)/3600, P_load(1:end-1)/1000, '-b', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('${P}_{l}$ (kW)', 'fontsize',15,'Interpreter','latex')
%title('Secondary tasks power', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 
subplot(3,2,6)
plot(t/3600, P_heat/1000, '-b', 'LineWidth', 2)
%plot(t/3600, G2_th*(T_bat - T_fluid)/1000, '-r', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$P_h$ (kW)', 'fontsize',15,'Interpreter','latex')
%title('Heater power', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 

figure()
subplot(2,1,1)
plot(t(1:end-1)/3600, -P_b(1:end-1)/1000 + U_bat^2/(2*R_bat)/1000*(1 -sqrt(1- 4*R_bat/U_bat^2*P_out(1:end-1))) , '-b', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$P_{out}$ - $P_b$ check (kW)', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on 
subplot(2,1,2)
plot(t(1:end-1)/3600, - P_m(1:end-1)/1000 + max(P_b(1:end-1)/1000, eff*P_b(1:end-1)/1000), '-b', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$P_m$ - $P_{b}$ check (kW)', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
grid on