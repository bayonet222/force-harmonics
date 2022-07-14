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

machine.s = 384;                                % Number of slots
machine.p = 160;                                % Number of pole pairs
machine.Nseg = 16;                              % Number of segments
machine.N_ph = 3;                               % Number of phases
machine.q = machine.s/machine.p/machine.N_ph/2; % Slots/pole/phase
machine.pf = 0.85;                              % Power factor

% Geometric parameters
machine.R_s = 4.89;                             % Stator inner radius [m]
machine.R_si = machine.R_s;
machine.R_so = 5.0;                             % Stator outer radius [m]
machine.R_r = 4.85;                             % Rotor radius [m];
machine.L = 1.7;                                % Machine length [m]

machine.g = 0.01;                               % Air gap [m]

machine.theta_s = 2*pi / machine.s;             % Slot pitch

machine.w_gap = 0.01;                           % Segment gap width UPDATE
machine.theta_gap = machine.w_gap/machine.R_s;  % Segment gap angle
machine.gap_loc = 19;                           % Tooth # of gap, from a-axis

machine.w_bi = 0.0466;                          % Back iron thickness [m]
machine.w_th = 0.0388;                          % Tooth width [m]
machine.theta_o = machine.theta_s - ...
    machine.w_th/machine.R_si;                  % Slot opening

% Magnetic parameters
machine.l_m = 0.03;                             % Magnet thickness [m]
machine.Br = 1.2;                               % Magnet remanence [T]
machine.alpha_m = 0.7;                          % Magnet span
machine.mu_r = 1.05;

% Winding parameters
machine.V_LL = 100e3 / 1.35;                    % Line-to-line voltage
machine.V_ph = machine.V_LL/sqrt(3);            % Phase voltage
machine.i_RMS = -machine.Pe/machine.V_ph/3/machine.pf;  % Phase current

machine.N_c = machine.s/2;                      % Number of coils
machine.n_s = 61;                               % Number of turns per coil
machine.I_s = machine.n_s * machine.i_RMS;      % Slot current

% Power Electronic parameters
machine.I_h_pu = [1 0.02 0.02];                      % Current time harmonics in pu
machine.h_PE = [1 37 41];                          % Time harmonic orders in the current

% Mechanical parameters
machine.rho_c = 7700;                           % Core density [kg/m3]
machine.k_i = 0.96;                             % Stacking factor

machine.E_plane = 180e9;                        % Young's modulus (in plane) [Pa]
machine.E_z = 100e9;                            % Young's modulus (axial) [Pa]
machine.nu_plane = 0.2;                        % Poisson ratio (in plane)
machine.nu_z = 0.1;                            % Poisson ratio (out of plane)

% Frame parameters
machine.h_f = 0.05;                             % Frame thickness [m]
machine.L_f = 2.0;                              % Frame length [m]
machine.rho_f = 7850;                           % Frame density [kg/m3] STEEL ASSUMED
machine.E_f = 200e9;                            % Frame Young's modulus [Pa]
machine.nu_f = 0.30;                            % Frame Poisson ratio

% Winding parameters
machine.rho_w = 8890;                           % Winding density

% Derived geometric parameters
machine.R_m = machine.R_r + machine.l_m;        % Magnet outer radius [m]
machine.d_s = machine.R_so - machine.R_s - machine.w_bi; % Slot depth

machine.D_f = 2*(machine.R_so + machine.h_f);   % Frame outer diameter [m]

% Derived electrical parameters
machine.f_e = machine.n_rpm/60 * machine.p;     % Electrical frequency [Hz]

% Calculate mass of the machine
[machine.M_c, machine.M_t, machine.M_w, machine.M_f] = masses(machine);

save params machine mu_0