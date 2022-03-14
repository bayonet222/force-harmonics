% Main program to analyse the harmonics in the ModHVDC machine
% Casper Klop 2022

% ------------------------------------------------------------------
%                      Load machine parameters
% ------------------------------------------------------------------
load params


% ------------------------------------------------------------------
%                    Define additional variables
% ------------------------------------------------------------------
parts = machine.Nseg;                               % Section of the machine analysed

N_theta_teeth = 101;                                % Number of samples/tooth
N_theta = (N_theta_teeth-1)*machine.s/parts + 1;    % Number of samples/segment
N_lambda = 25;                                      % Number of Fourier components

theta_vect = linspace(0, 2*pi/parts, N_theta);      % Theta [rad]
theta_vect_teeth = linspace(0, 2*pi/machine.s, N_theta_teeth);    

theta_r = linspace(0, 2*pi/parts, 1000);            % Rotor position
t_vect = theta_r/machine.omega_m;                   % Time [s]

% ------------------------------------------------------------------
%                      Calculate magnetic field
% ------------------------------------------------------------------

% Magnetic field due to Permanent Magnets
[B_mmu, B_PM_sl] = B_PM_slotless(theta_vect, t_vect, machine);

% Magnetic field due to armature reaction
[B_mv, B_arm_sl] = B_arm_slotless(theta_vect, t_vect, machine);

% Total slotless magnetic field, superposition of PM and armature
% B_slotless = B_PM_sl + B_arm_sl;
B_slotless = B_PM_sl + 0.0 * B_arm_sl;           % Currently focusing on B_PM at as armature reaction too strong

% Use Zarko's method for relative permeance along 1 slot or gap
[lambda_slot, lambda_an_slot, lambda_bn_slot] = Zarko(machine, theta_vect_teeth, N_lambda, 'slot');
[lambda_gap, lambda_an_gap, lambda_bn_gap] = Zarko(machine, theta_vect, N_lambda*machine.s/machine.Nseg, 'gap');

% Repeat the slot permeance s/Nseg times and add to the gap permeance
lambda_slot_total = [repmat(lambda_slot(1:end-1), 1, machine.s/machine.Nseg), lambda_slot(end)];
lambda = lambda_slot_total + lambda_gap - 1;

% Obtain the slotted flux density
B = B_slotless .* conj(lambda);

% Perform Fourier transform of flux densities
% TBD

% ------------------------------------------------------------------
%                          Force calculations
% ------------------------------------------------------------------

% Total complex forces and its FFT
f = B.^2 / (2 * mu_0);
[r, f_rr] = fft_ss(real(f(1,:)), 300, parts);

% Constant forces and torque
f_avg = transpose(mean(f, 2));
T_e = machine.R_s^2 * 2*pi * machine.L * imag(f_avg);

% Calculate zeroth radial order. Correct?
f_rr(1) = (max(real(f_avg)) - min(real(f_avg)))/2;

% Calculate cogging torque due to slots and segments
T_c = cogging_torque(B_PM_sl(1,:), lambda_slot_total, machine, parts, theta_r, 'slot');
T_c_seg = cogging_torque(B_PM_sl(1,:), lambda_gap, machine, parts, theta_r, 'gap');

% ------------------------------------------------------------------
%                        Mechanical calculations
% ------------------------------------------------------------------

% Calculate static deformations for modes r
Y_ms = deformation(machine, f_rr(1:20), r(1:20));

f_n = Natural_frequencies(machine, r(1:20));

% ------------------------------------------------------------------
%                          Plot results
% ------------------------------------------------------------------

% Defaults
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',12);
set(0, 'DefaultLineLineWidth', 2);

set(gcf,'color','w');
set(0,'DefaultFigureColormap',parula);

theta_ticks = 0:2*pi/parts/4:2*pi/parts;
theta_labels = strrep(string(sym(theta_ticks)), 'pi', '\pi');

% Show accuracy of Fourier series of permeance
[plt.lambda_re_slot, plt.lambda_im_slot] = ...
    Zarko_reconstruct(lambda_an_slot, lambda_bn_slot, machine.s, theta_vect_teeth, lambda_slot);
[plt.lambda_re_gap, plt.lambda_im_gap] = ...
    Zarko_reconstruct(lambda_an_gap, lambda_bn_gap, machine.Nseg, theta_vect, lambda_gap);

% Plot slotless fields
plt.B_PM_sl = figure;
plot(theta_vect, B_PM_sl(1,:));
title('Slotless radial flux density due to PMs')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_arm_sl_r = figure;
plot(theta_vect, real(B_arm_sl(1,:)));
title('Slotless radial flux density due to armature reaction')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_arm_sl_t = figure;
plot(theta_vect, imag(B_arm_sl(1,:)));
title('Slotless tangential flux density due to armature reaction')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot permeance 
plt.lambda_re = figure;
plot(theta_vect, real(lambda));
title('Real part of the relative permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.lambda_im = figure;
plot(theta_vect, imag(lambda));
title('Imaginary part of the relative permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot slotted flux density
plt.B_r = figure;
plot(theta_vect, real(B(1,:)));
title('Slotted radial flux density at t=0')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_t = figure;
plot(theta_vect, imag(B(1,:)));
title('Slotted tangential flux density at t=0')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot forces
plt.f_r = figure;
plot(theta_vect, real(f(1,:)));
title('Slotted radial forces at t=0')
xlabel('Angular position [rad]')
ylabel('Force density [N/m^2]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.f_r_fft = figure;
bar(r, f_rr, 'FaceAlpha', 1);
title('Harmonics of the radial force')
xlabel('Spatial order')
ylabel('Force density [N/m^2]')

plt.f_t = figure;
plot(theta_vect, imag(f(1,:)));
title('Slotted tangential forces at t=0')
xlabel('Angular position [rad]')
ylabel('Force density [N/m^2]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot cogging torque for slots and slots+segments
plt.Tc = figure;
plot(theta_r, T_c, theta_r, (T_c + T_c_seg));
title('Cogging torque')
xlabel('Rotor angular position [rad]')
ylabel('Torque [Nm]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot average forces
plt.f_avg_r = figure;
plot(t_vect, real(f_avg));
title('Average radial force')
xlabel('Time [s]')
ylabel('Force density [N/m^2]')
grid on;

plt.Te = figure;
plot(t_vect, T_e);
title('Total torque')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on;

% Plot deformation
plt.Yms = figure;
bar(r(1:length(Y_ms)), Y_ms, 'FaceAlpha', 1);
title('Static deformation for different modes')
xlabel('Spatial order')
ylabel('Deflection [m]')

% ------------------------------------------------------------------
%                            Save figures
% ------------------------------------------------------------------

export_fig(plt.lambda_re_slot, 'Figures/lambda_re_slot', '-eps', '-transparent')
export_fig(plt.lambda_im_slot, 'Figures/lambda_im_slot', '-eps', '-transparent')
export_fig(plt.lambda_re_gap, 'Figures/lambda_re_gap', '-eps', '-transparent')
export_fig(plt.lambda_im_gap, 'Figures/lambda_im_gap', '-eps', '-transparent')
export_fig(plt.B_PM_sl, 'Figures/B_PM_sl', '-eps', '-transparent')
export_fig(plt.B_arm_sl_r, 'Figures/B_arm_sl_r', '-eps', '-transparent')
export_fig(plt.B_arm_sl_t, 'Figures/B_arm_sl_t', '-eps', '-transparent')
export_fig(plt.lambda_re, 'Figures/lambda_re', '-eps', '-transparent')
export_fig(plt.lambda_im, 'Figures/lambda_im', '-eps', '-transparent')
export_fig(plt.B_r, 'Figures/B_r', '-eps', '-transparent')
export_fig(plt.B_t, 'Figures/B_t', '-eps', '-transparent')
export_fig(plt.f_r, 'Figures/f_r', '-eps', '-transparent')
export_fig(plt.f_r_fft, 'Figures/f_r_ftt', '-eps', '-dNOSAFER')
export_fig(plt.f_t, 'Figures/f_t', '-eps', '-transparent')
export_fig(plt.Tc, 'Figures/Tc', '-eps', '-transparent')
export_fig(plt.f_avg_r, 'Figures/f_avg_r', '-eps', '-transparent')
export_fig(plt.Te, 'Figures/Te', '-eps', '-transparent')
export_fig(plt.Yms, 'Figures/Yms', '-eps')

% ------------------------------------------------------------------
%                             Create GIFs
% ------------------------------------------------------------------

%make_gif(theta_vect, t_vect, real(f), 'ForceRadial')