function T = getTemperature(t, nu, t_0, thermal_model, cos_phi_s)
    
    Q = getHeat(nu, thermal_model, cos_phi_s); % compute solar heat
    T=zeros(size(t)); 
    T(1)=t_0;
    for i=1:length(t)-1 
    %T(i+1)=T(i) + thermal_model.d/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*T(i).^4);
    %f=1/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*T(i).^4);

    k1 = 1/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*T(i).^4);
    k2 = 1/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*(T(i)+k1/2*thermal_model.d).^4);
    k3 = 1/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*(T(i)+k2/2*thermal_model.d).^4);
    k4 = 1/thermal_model.C_th*(Q(i) - thermal_model.epsilon*thermal_model.sigma*thermal_model.S*(T(i)+k3*thermal_model.d).^4);
    T(i+1)=T(i) + 1/6*thermal_model.d*(k1 + 2*k2 + 2*k3 + k4); 
    end 
end

