function Y_ms = deformation(machine_params, f_rm, m)
% Returns the static deformations for forces at mode m
%   Based on the geometric parameters of the machine, the 
%   amplitude of the forces and the mode shapes, this function
%   returns the static deformations. Based on Cassoret 2011
%   The first force magnitude should represent the 0th order
R = machine_params.R_s;                             % Stator inner radius
R_y = machine_params.R_so - machine_params.w_bi/2;  % Average yoke radius
E = machine_params.E_plane;                         % Young's modulus
T_y = machine_params.w_bi;                          % Back iron width

Y_0s = R * R_y * f_rm(1) / (E * T_y);               % Zeroth order deformation
Y_ms = [Y_0s ...
        12 * R * R_y^3 * f_rm(2:end) ./ ...
        (E * T_y^3 * (m(2:end).^2 - 1).^2)];
end

