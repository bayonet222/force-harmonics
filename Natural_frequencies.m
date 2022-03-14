function [f_mc] = Natural_frequencies(machine_params, m)
%Estimate the natural frequencies for the stator
%   Gives the natural frequencies for the stator core and
%   the stator-frame system for modes r. Based on Gieras Ch5

s = machine_params.s;                           % Number of slots
h_c = machine_params.w_bi;                      % Stator thickness, assumed to be the yoke thickness
D_c = 2*machine_params.R_so - h_c;              % Mean diameter                    
rho_c = machine_params.rho_c;                   % Core density
k_i = machine_params.k_i;                       % Stacking factor
L_i = machine_params.L;                         % Length
h_sl = machine_params.d_s;                      % Slot depth
w_th = machine_params.w_th;                     % Tooth width
E_c = machine_params.E_plane;                   % Young's modulus
nu_c = machine_params.nu_plane;                 % Poisson ratio

kappa2 = h_c^2 / (3 * D_c^2);                   % Nondimensional thickness


% Calculate masses of different components
M_c = pi*D_c*h_c*L_i * rho_c * k_i;             % Stator core mass 
M_t = s * h_sl * w_th * L_i * rho_c * k_i;      % Teeth mass

% Calculate mass addition factors
k_md = 1;                                       % For displacement

M_0 = M_c * k_md;


% Calculate stiffnesses
% Eq. (5.25)
K_mc = 4 * roots(kappa2, m).^2 / D_c * pi * L_i * h_c * E_c / (1-nu_c^2);
M_m = M_c;

f_mc = 1/(2*pi) * sqrt(K_mc ./ (M_m + M_t));

end

function Omega_m = roots(kappa2, m)
    % Roots of the second order characteristic equation
    
Omega_m = 0.5 * sqrt((1 + m.^2 + kappa2*m.^4) + ...
    sqrt((1 + m.^2 + kappa2 * m.^4).^2 - 4 * kappa2 * m.^6));
Omega_m(1) = 1;
end
