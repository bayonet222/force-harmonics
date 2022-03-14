% ---------------------------------------------------------------------
%                              Constants
% ---------------------------------------------------------------------
global mu_0 
mu_0 = 4 * pi * 10^-7;                          % Vacuum permeability

% ---------------------------------------------------------------------
%                          Machine parameters
% ---------------------------------------------------------------------
machine.Pe = 10e6;                              % Power [W]

machine.n_rpm = 10;                             % Rotational speed [RPM]
machine.omega_m = machine.n_rpm*2*pi/60;        % Rotational speed [rad/s]

machine.T = machine.Pe/machine.omega_m;         % Torque [Nm]

machine.s = 192;                                % Number of slots
machine.p = 80;                                 % Number of pole pairs
machine.Nseg = 16;                              % Number of segments
machine.N_ph = 3;                               % Number of phases
machine.q = machine.s/machine.p/machine.N_ph/2; % Slots/pole/phase

% Geometric parameters
machine.R_s = 3.517;                            % Stator inner radius [m]
machine.R_si = machine.R_s;
machine.R_so = 3.7;                             % Stator outer radius [m]
machine.R_r = 3.477;                            % Rotor radius [m];
machine.L = 1.5;                                % Machine length [m]

machine.g = 0.01;                               % Air gap [m]

machine.theta_s = 2*pi / machine.s;             % Slot pitch
machine.theta_o = machine.theta_s/2;            % Slot opening UPDATE

machine.w_gap = 0.01;                           % Segment gap width UPDATE
machine.theta_gap = machine.w_gap/machine.R_s;  % Segment gap angle

machine.w_bi = 0.069;                           % Back iron thickness [m]
machine.w_th = 0.058;                           % Tooth width [m]

% Magnetic parameters
machine.l_m = 0.03;                             % Magnet thickness [m]
machine.Br = 1.2;                               % Magnet remanence [T]
machine.alpha_m = 0.7;                          % Magnet span
machine.mu_r = 1.05;

% Winding parameters
machine.V_LL = 100e3 / 1.35;                    % Line-to-line voltage
machine.V_ph = machine.V_LL/sqrt(3);            % Phase voltage
machine.i_RMS = machine.Pe/machine.V_ph/3;      % Phase current

machine.N_c = machine.s/2;                      % Number of coils
machine.n_s = 187;                              % Number of turns per coil
machine.I_s = machine.n_s * machine.i_RMS;      % Slot current

% Mechanical parameters
machine.rho_c = 7700;                           % Core density [kg/m3]
machine.k_i = 0.96;                             % Stacking factor

machine.E_plane = 180e9;                        % Young's modulus (in plane) [Pa]
machine.E_x = 100e9;                            % Young's modulus (axial) [Pa]
machine.nu_plane = 0.22;                        % Poisson ratio (in plane)
machine.nu_x = 0.64;                            % Poisson ratio (out of plane)



% Derived geometric parameters
machine.R_m = machine.R_r + machine.l_m;        % Magnet outer radius [m]
machine.d_s = machine.R_so - machine.R_s - machine.w_bi; % Slot depth

save params machine mu_0