function B_arm = B_arm_slotless_PE(theta_vect, t_vect, machine_params)
%Function B_arm of theta, taking time harmonics into account
%   Returns the function of the magnetic flux density
%   due to the armature reaction at the stator teeth
%   as function of theta

global mu_0

% Load machine parameters
p = machine_params.p;                               % Number of pole pairs
s = machine_params.s;                               % Number of slots
theta_o = machine_params.theta_o;                   % Slot opening pitch
N_t = machine_params.n_s * machine_params.N_c ...
    /machine_params.N_ph;                           % Turns per phase
R_s = machine_params.R_s;                           % Inner stator radius
R_r = machine_params.R_r;                           % Rotor radius
I_1 = machine_params.i_RMS * sqrt(2);               % Amplitude of fundamental current
omega_m = machine_params.omega_m;                   % Rotation speed

% Calculate the current magnitudes (I_vect given in pu)
I_h = I_1 * machine_params.I_h_pu;
h_PE = machine_params.h_PE;

% Spatial orders, odd, non-triplen multiples of F1
F_1 = gcd(p, s);                                    % Number of parts
v_seg = 1:2:50;
v_seg(~mod(v_seg, 3)) = [];

v = transpose(v_seg) * F_1;

% Create theta-t matrix and loop through time
B_arm = zeros(length(t_vect), length(theta_vect));

for k = 1:length(h_PE)
% Loop through current harmonics

    % Radial and tangential magnitudes
    B_mv_r = 3*mu_0/2 * J_v(N_t, R_s, v, s, theta_o) .* F_v(R_r, R_s, v) * I_h(k);
    B_mv_t = 3*mu_0/2 * J_v(N_t, R_s, v, s, theta_o) .* G_v(R_r, R_s, v) * I_h(k);

    for i = 1:length(t_vect)
        B_arm(i,:) = B_arm(i,:) + sum(B_mv_r .* sin(rot(v, F_1).* v.*theta_vect - h_PE(k) * p*omega_m*t_vect(i)) ...
        + 1j * B_mv_t .* rot(v, F_1) .* cos(rot(v, F_1).* v.*theta_vect - h_PE(k) * p*omega_m*t_vect(i)));
    end
    
end

end

function J_func = J_v(N_t, R_s, v, s, theta_o)
% J_m according to Zhu (2010)

J_func = 2/pi * N_t/R_s * K_sov(v, theta_o) .* K_pv(v, s) .* K_dv(v);
end

function slot_opening = K_sov(v, theta_o)
% The slot opening factor from Zhu (2010)

slot_opening = sin(v * theta_o/2) ./ (v * theta_o/2);
end

function pitch_factor = K_pv(v, s)
% The pitch factor as described by Zhu and Gieras
% K_pv = sin(v * beta_s/2), beta in mechanical radians
pitch_factor = sin(v * pi /s);

end

function distribution_factor = K_dv(v)
% The distribution factor is equal to 1
distribution_factor = 1 .^v;
end

function rotation_sign = rot(v, F1)
% The rotation factor according to Yokoi (2016)

k_taum = round(sin(pi * (v/F1-1)/2),5) ./ round((3 * sin(pi*(v/F1-1)/6)),5);
k_taum(isnan(k_taum))=1;
k_taup = round(sin(pi * (v/F1+1)/2),5) ./ round((3 * sin(pi*(v/F1+1)/6)),5);
k_taup(isnan(k_taup))=1;

rotation_sign = sign(k_taup-k_taum);
end

function geom = F_v(R_r, R_s, v)
% The geometrical function of Zhu (2010) for r = R_s
% should approximate 1/v * Rs/g

geom = (1 + (R_r/R_s).^(2*v)) ./ (1 - (R_r/R_s).^(2*v));
end

function geom2 = G_v(R_r, R_s, v)
% Tangential geometric function of Zhu (2010) equals 1 
% for r = R_s

geom2 = ones(length(v), 1);
end