% parameters
beta_s = 2*pi / 192;
beta_o = beta_s/2;
b14 = 3.517 * beta_o;
t1 = beta_s * 3.517;
g = 0.01;
h_m = 0.03;
muR = 1.05;

v = [1 5 7 11 13 17 19 23 25 29]*16;     % Spatial orders

% Slot opening factors
k_so_zhu = sin(v * beta_o/2) ./ (v * beta_o/2);

kappa = b14 / (g + h_m/muR);
rho = kappa/(5+kappa) * 2 * sqrt(1+kappa^2) / (sqrt(1+kappa^2) -1);
k_so_gieras = sin(v * rho * pi * b14/2/t1) ./ (v * rho * pi * b14/2/t1);

% Slot pitch factors
k_p_zhu = sin(v * beta_s/2);

gamma_s = pi / 0.4/3;
ep = pi-gamma_s;
k_p_me = cos(v*ep/2);

% Distribution factors
k_d_zhu = sin(v * pi*80/192) .* sin(v/16 * pi/2);

k_d_me = ones(1, length(v));

% Zhu factor
Fm = (1 + (3.477/3.517).^(2*v/16)) ./ (1 - (3.477/3.517).^(2*v/16));
% This approx equals R_s / g * 1/v. Harmonic order 

figure
plot(v, Fm)