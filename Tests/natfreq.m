m=2:20;
E_c = 180E9;
rho_c = 7700;
L = 1.7;
h_c = 0.0466;
D_c = 10-h_c;
nu_c = 0.22

I_c = h_c^3 * L/12;

% First method, Hoppe
f_1 = 2/(pi*D_c^2) * m.*(m.^2-1)./sqrt(m.^2 +1) * sqrt(E_c*I_c/rho_c/L/h_c)

% Second method, DM
kappa2 = h_c^2/(3*D_c^2);

% TEST MINUS SIGN HERE
Omega_m = 0.5*sqrt((1+m.^2 + kappa2*m.^4) - ...
    sqrt((1 + m.^2 + kappa2*m.^4).^2 - 4*kappa2*m.^6));

f_2 = Omega_m/pi/D_c * sqrt(E_c/(rho_c * (1-nu_c^2)))