function [I_ds, I_L, V_t]   = getSolarArrayParameters(T, d, nu, param)
    
    % SA cell parameters 
    V_oc = temp(param.V_oc_cell_SA, param.T_ground_SA, param.alpha_V_oc_cell_SA, T); % open circuit voltage cell (A)
    I_sc = temp(param.I_sc_cell_SA, param.T_ground_SA, param.alpha_I_sc_cell_SA, T); % short circuit current cell (A)  
    V_t = param.k_B*T/ param.q_e; % thermal voltage cell (V) 
    
    % Currents: 
    I_ds = param.N_p_cell_SA * I_sc .* exp(-V_oc./(V_t)); % diode saturation current for the whole SA (A)
    I_L = param.N_p_cell_SA*  param.I_0*param.G_0*(param.R_sun/d)^2*nu; % irradiation current 
    V_t = V_t*param.N_s_cell_SA; % thermal voltage of the array 

end

