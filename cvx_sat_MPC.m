tic 
%% Initialise 
T_bat_MPC = zeros(L,1); 
P_out_MPC = zeros(L,1);
P_load_MPC = zeros(L,1); 
P_heat_MPC = zeros(L,1);
P_m_MPC = zeros(L,1);
E_MPC = zeros(L,1);
P_b_MPC = zeros(L,1);

E_max = 830822400/150*10*2*2;
E_min = 0.3*E_max; 
E_0 = E_min;
T_0 = T_0_bat;
avg=0;

for i=1:L
i
%% Resize problem
N = min(L-i+1, L_local); 
Phi=ones(N,1);
a_bat=(1 - exp(-G_th*t(i:i+N-1)/C_th_bat))/G_th ;
b_bat = T_0*exp(-G_th*t(i:i+N-1)/C_th_bat);


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
           P_out*P + V_d(i:i+N-1,1).*I_s(i:i+N-1,1)  ==  P_load*P + P_a(i:i+N-1,1) + P_HET(i:i+N-1,1)  +  P_heat*P + P_comp/2  %I_L - I_ds*(exp(V_x)-1.0) + abs(P_heat*P)/U_heat
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
           a_bat.*(P_comp+P_heat*P + G1_th*T_surface(i:i+N-1,1) + G2_th*T_fluid) + b_bat <= Phi*(T_max_bat)
           T_min_bat*Phi <= a_bat.*( P_comp+ P_heat*P + G1_th*T_surface(i:i+N-1,1) + G2_th*T_fluid) + b_bat
 
cvx_end
cvx_status

% T_bat = a_bat.*(P_comp+P_heat*P + G1_th*T_surface(i:i+N-1,1) + G2_th*T_fluid) + b_bat;
% P_out = P_out*P; P_load = P_load*P;  P_heat = P_heat*P;   %I_heat = I_heat*I; I_out = I_out*I;
% P_m=P_m*P; E= E_*E; P_b = P_b*P;
% t = t(i:i+N-1,1);
% P_s = V_d(i:i+N-1,1).*I_s(i:i+N-1,1);

%% Rescaling and store first element
T_bat = a_bat.*(P_comp+P_heat*P + G1_th*T_surface(i:i+N-1,1) + G2_th*T_fluid) + b_bat;
T_bat_MPC(i) = T_bat(1); 
P_out_MPC(i) = P_out(1)*P; P_load_MPC(i) = P_load(1)*P;  P_heat_MPC(i) = P_heat(1)*P;   %I_heat = I_heat*I; I_out = I_out*I;
P_m_MPC(i)=P_m(1)*P; E_MPC(i)= E_*E(1); P_b_MPC(i) = P_b(1)*P;


%% MPC update:
if(i < L)
    if (P_b(1) < 0)
        E_0 = E_*E(1) -d*P_b(1)*P*eff/0.7;
    else 
        E_0 = E_*E(1) -d*P_b(1)*P;
    end
% if (E_0 < E_min)
%    E_0 = E_min;  
% end
E_0 = E(1)*E_ -d*max(P_b(1)*P, P_b(1)*P*eff); %eff 
T_0 = T_bat(2); 
end

%% Test with model uncertainty
if(i==1)
    Energy = E_*E;
    Power = P_b*P; 
    heat = P_heat*P; 
end 
if (i ==  L_local*3) 
    break 
end 
end 
toc

%% Rename variables 
T_bat = T_bat_MPC(1:i);
P_out = P_out_MPC(1:i); P_load = P_load_MPC(1:i);  P_heat = P_heat_MPC(1:i); 
P_m=P_m_MPC(1:i); E= E_MPC(1:i); P_b = P_b_MPC(1:i);
P_s = V_d(1:i).*(I_L(1:i) - I_ds(1:i).*(exp(V_x(1:i))-1.0)); 
t = t(1:i); 

 

%% Plot 
% figure()
% subplot(2,2,1)
% plot(t/3600, P_s, '-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$P_s$ (A)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% subplot(2,2,3)
% plot(t/3600, T_bat-273.15, '-b', 'LineWidth', 2)
% hold on 
% plot(t/3600, ones(size(t))*(T_max_bat-273.15), '--r', 'LineWidth', 2)
% plot(t/3600, ones(size(t))*(T_min_bat-273.15), '--r', 'LineWidth', 2)
% axis([0 t(end)/3600 (T_min_bat -273.15)*0.8 (T_max_bat-273.15)*1.2])
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$T_b$ ($^\circ$C)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% subplot(2,2,2)
% %plot(t/3600, T_surface -273.15, '-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$T_c$ ($^\circ$C)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% subplot(2,2,4)
% hold on
% plot(t/3600, P_heat/1000/5, '-b', 'LineWidth', 2)
% plot(t/3600, G2_th*(T_bat - T_fluid)/1000/5, '-r', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$P_h$ (kW)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% 
% 
% figure()
% subplot(2,2,1)
% plot(t(1:end-1)/3600,  E(1:end-1) -E(2:end) -d*P_m(1:end-1), '-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('SOC check (\%)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% subplot(2,2,3)
% plot(t(1:end-1)/3600, - P_m(1:end-1) + max(P_b(1:end-1), eff*P_b(1:end-1)), '-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$P_m$ - $P_{b}$ check (A)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on
% subplot(2,2,2)
% plot(t/3600,P_out + V_d.*I_s - P_load -P_a -P_HET -  P_heat -3500,'-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('Power balance check (W)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 
% subplot(2,2,4)
% plot(t/3600, -P_b + U_bat^2/(2*R_bat)*(1 -sqrt(1- 4*R_bat/U_bat^2*P_out)) , '-b', 'LineWidth', 2)
% xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
% ylabel('$P_out$ - $P_b$ check (W)', 'fontsize',15,'Interpreter','latex')
% set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
% grid on 


%% Final figures: 
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

figure()
plot(t/3600, E/1e+6, '-b', 'LineWidth', 2)
hold on 
for i=1:L_local-1
    if (Power(i) < 0)
        E(i+1) = E(i) -d*Power(i)*eff/0.7;
    else 
    E(i+1) = E(i) -d*Power(i);
    end
end
plot(t/3600, E/1e+6, '--k', 'LineWidth', 2)
plot(t/3600, ones(size(t))*E_max/1e+6, '--r', 'LineWidth', 2)
plot(t/3600, ones(size(t))*E_min/1e+6, '--r', 'LineWidth', 2)
xlabel('Time (h)', 'fontsize',15,'Interpreter','latex')
ylabel('$E$ (MJ)', 'fontsize',15,'Interpreter','latex')
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)
legend({'MPC', 'Single-shot'}, 'fontsize',15,'Interpreter','latex')
grid on

