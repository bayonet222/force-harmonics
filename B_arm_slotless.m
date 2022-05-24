function B_arm = B_arm_slotless(theta_vect, t_vect, machine_params)
%Function B_arm of theta
%   Returns the function of the magnetic flux density
%   due to the armature reaction at the stator teeth
%   as function of theta

global mu_0

p = machine_params.p;
s = machine_params.s;
theta_o = machine_params.theta_o;
N_t = machine_params.n_s * machine_params.N_c ...
    /machine_params.N_ph;                           % Turns per phase
R_s = machine_params.R_s;
R_r = machine_params.R_r;
I_1 = machine_params.i_RMS * sqrt(2);               % Amplitude of fundamental current

omega_m = machine_params.omega_m;                   % Rotation speed

F_1 = gcd(p, s);                                    % Number of parts
v = transpose([1 5 7 11 13 17 19 23 25 29 31 35 37 41 43 47 49]*F_1);    % Spatial orders

B_mv_r = 3*mu_0/2 * J_v(N_t, R_s, v, s, theta_o) .* F_v(R_r, R_s, v) * I_1;
B_mv_t = 3*mu_0/2 * J_v(N_t, R_s, v, s, theta_o) .* G_v(R_r, R_s, v) * I_1;

% Create theta-t matrix and loop through time
B_arm = zeros(length(t_vect), length(theta_vect));

for i = 1:length(t_vect)
    B_arm(i,:) = sum(B_mv_r .* sin(rot(v, F_1).* v.*theta_vect - p*omega_m*t_vect(i)) ...
    + 1j * B_mv_t .* rot(v, F_1) .* cos(rot(v, F_1).* v.*theta_vect - p*omega_m*t_vect(i)));

%     B_arm(i,:) = sum(B_mv_r .* sin(v.*theta_vect - p*omega_m*t_vect(i)) ...
%     + 1j * B_mv_t .* cos(v.*theta_vect - p*omega_m*t_vect(i)));
end

end

function J_func = J_v(N_t, R_s, v, s, theta_o)
% J_m according to Zhu (2010)

J_func = 2/pi * N_t/R_s * K_sov(v, theta_o) .* K_pv(v, s) .* K_dv(v);
end

function slot_opening = K_sov(v, theta_o)
% The slot opening factor from Zhu2010

slot_opening = sin(v * theta_o/2) ./ (v * theta_o/2);
end

function pitch_factor = K_pv(v, s)
% The pitch factor as described by Zhu and Gieras
% K_pv = sin(v * beta_s/2), beta in mechanical radians
pitch_factor = sin(v * pi /s);

end

function distribution_factor = K_dv(v)
% The distribution factor according to Zhu/
% In fact, it only creates a [1 -1 1 -1 ...] factor for v=[1 3 5 7]

%distribution_factor = sin(v/F1 * pi/2);
distribution_factor = 1 .^v;
end

function rotation_sign = rot(v, F1)

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