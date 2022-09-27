function Q = getHeat(nu, th, cos_phi_s)
    
    G = th.G_0*(th.R_sun/th.AU)^2.*max(0, cos_phi_s); % intensity of solar irradiation at surface
    
    Q = nu*th.gamma.*G*th.S*cos(th.phi)+237.0*th.S*th.epsilon*(th.R_earth/(th.z))^2; 

end

