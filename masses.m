function [M_c, M_t, M_w, M_f] = masses(machine_params)
%Calculate masses of the system

s = machine_params.s;                           % Number of slots
h_c = machine_params.w_bi;                      % Stator thickness, assumed to be the yoke thickness
D_c = 2*machine_params.R_so - h_c;              % Mean diameter                    
rho_c = machine_params.rho_c;                   % Core density
k_i = machine_params.k_i;                       % Stacking factor
L_i = machine_params.L;                         % Length
h_sl = machine_params.d_s;                      % Slot depth
w_th = machine_params.w_th;                     % Tooth width
w_s = 2*pi/s * machine_params.R_si - w_th;      % Slot width

% Frame parameters
D_f = machine_params.D_f;                       % Frame outer diameter
L_f = machine_params.L_f;                       % Frame length
h_f = machine_params.h_f;                       % Frame thickness
rho_f = machine_params.rho_f;                   % Frame density
R_f = 0.5 * (D_f - h_f);                        % Mean frame radius

% Winding parameters
rho_w = machine_params.rho_w;                   % Winding density

% Calculate masses of different components
M_c = pi*D_c*h_c*L_i * rho_c * k_i;             % Stator core mass 
M_t = s * h_sl * w_th * L_i * rho_c * k_i;      % Teeth mass
M_w = (s * w_s * h_sl * L_i + ...
    (w_s + w_th) * D_c * pi * h_sl)* rho_w;     % Winding mass (slot + overhang)
M_f = pi*2*R_f*h_f*L_f * rho_f;                 % Frame mass
