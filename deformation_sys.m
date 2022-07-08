function A_m = deformation_sys(machine_params, f_rm, freq_sys)
% Generalized expression for the deformation of the system 
%   Taking into account the frame and windings
R = machine_params.R_so;                        % Outer stator radius
L = machine_params.L_f;                         % Frame length

% Machine mass
M_c = machine_params.M_c;                       % Stator core mass
M_t = machine_params.M_t;                       % Teeth mass
M_w = machine_params.M_c;                       % Winding mass (slot + overhang)
M_f = machine_params.M_f;                       % Frame mass
M = M_c + M_t + M_w + M_f;                      % Total mass

A_m = R * L ./ (M * 2*pi * freq_sys.^2) .* f_rm;
end

