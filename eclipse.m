function nu = eclipse(state)
rx = state(1,1);
ry = state(2,1); 
R_sun =  696e+3 ; 
R_earth = 6378.0 ;
AU = 149597870 ; 
r0 = 42164; 
r_sun_earth = [-AU ; 0]; 
r_sat_earth = [rx; ry]; 
chi = R_earth / (R_sun -R_earth)* norm(r_sun_earth); 
a_e = asin(R_earth/chi);
theta = r_sun_earth'*r_sat_earth / (norm(r_sun_earth)*norm(r_sat_earth));

% a = pi - acos (theta) ; 
% gamma = (2*chi*cos(a_e) - sqrt( (2*chi*cos(a_e))^2 - 4*(chi^2-r0^2) )) /2; 
% a_e2 = acos(  (gamma^2 - r0^2 - chi^2) /(-2*r0*chi)); 
% condition : abs(a) < a_e2 

s = [1; 0]; 
r_proj_earth = (s'*r_sat_earth)*s; 
delta = r_sat_earth - r_proj_earth;  
xi = (chi -  norm(r_proj_earth))*tan(a_e);

%% Umbra 
if (norm(delta) <= xi && (s'*r_sat_earth) > 0)
    nu_u = 0;  
else 
    nu_u = 1; 
end 

%% Penumbra 
chi_p = R_earth / (R_sun + R_earth)* norm(r_sun_earth); 
a_p = asin(R_earth/chi_p);
kappa = (chi_p + norm(r_proj_earth))*tan(a_p);

if (norm(delta) <= kappa && (s'*r_sat_earth) > 0)
    nu_p = (norm(delta)  -xi)/(kappa-xi);  
else 
    nu_p = 1; 
end 

nu = nu_u*nu_p;

end 

