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

% Stationary part
N_theta_teeth = 101;                                % Number of samples/tooth
N_theta = (N_theta_teeth-1)*machine.s/parts + 1;    % Number of samples/segment
N_lambda = 25;                                      % Number of Fourier components

theta_vect = linspace(0, 2*pi/parts, N_theta);      % Theta [rad]
theta_vect_teeth = linspace(0, 2*pi/machine.s, N_theta_teeth);    

% Gap index
gap_idx = (machine.gap_loc - machine.s/parts/2) * (N_theta_teeth-1);

% Rotating part
theta_r = linspace(0, 2*pi/parts, 1001);            % Rotor position
t_vect = theta_r/machine.omega_m;                   % Time [s]

% ------------------------------------------------------------------
%                      Calculate magnetic field
% ------------------------------------------------------------------

% Magnetic field due to Permanent Magnets
B_PM_sl = B_PM_slotless(theta_vect, t_vect, machine);

% Magnetic field due to armature reaction
B_arm_sl = B_arm_slotless(theta_vect, t_vect, machine);

% Total slotless magnetic field, superposition of PM and armature
B_slotless = B_PM_sl + B_arm_sl;          

% Use Zarko's method for relative permeance along 1 slot or gap
[lambda_slot, lambda_an_slot, lambda_bn_slot] = Zarko(machine, theta_vect_teeth, N_lambda, 'slot');
[lambda_gap_c, lambda_an_gap, lambda_bn_gap] = Zarko(machine, theta_vect, N_lambda*machine.s/machine.Nseg, 'gap');

% Shift the segment gap to the correct location
lambda_gap = circshift(lambda_gap_c, gap_idx);

% Repeat the slot permeance s/Nseg times and add to the gap permeance
lambda_slot_total = [repmat(lambda_slot(1:end-1), 1, machine.s/machine.Nseg), lambda_slot(end)];
lambda = lambda_slot_total + lambda_gap - 1;

% Obtain the slotted flux density
B_noload = B_PM_sl .* conj(lambda);
B_arm = B_arm_sl .* conj(lambda);
B_unseg = B_slotless .* conj(lambda_slot_total);
B = B_slotless .* conj(lambda);

% Perform Fourier transform of flux densities
[k, b_rk_NL] = fft_ss(real(B_noload(1, :)), 100, parts);
[~, b_tk_NL] = fft_ss(imag(B_noload(1, :)), 100, parts);

[~, b_rk_arm] = fft_ss(real(B_arm(1, :)), 100, parts);
[~, b_tk_arm] = fft_ss(imag(B_arm(1, :)), 100, parts);

[~, b_rk] = fft_ss(real(B(1, :)), 100, parts);
[~, b_tk] = fft_ss(imag(B(1, :)), 100, parts);
[B_fx, B_ft, B_2d_fft] = fft_2D(real(B), parts, t_vect(end));

% ------------------------------------------------------------------
%                          Force calculations
% ------------------------------------------------------------------

% Total complex forces and its FFT
f_noload = B_noload.^2 / (2 * mu_0);

f_unseg = B_unseg.^2 / (2 * mu_0);
[r, f_rr_unseg] = fft_ss(real(f_unseg(1,:)), 300, parts);

f = B.^2 / (2 * mu_0);
[~, f_rr] = fft_ss(real(f(1,:)), 300, parts);
[~, f_rt] = fft_ss(imag(f(1,:)), 300, parts);

[f_fx, f_ft, f_2d_fft] = fft_2D(real(f), parts, t_vect(end));

% Constant forces and torque
f_NL_avg = transpose(mean(f_noload, 2));
T_e_NL = machine.R_s^2 * 2*pi * machine.L * imag(f_NL_avg);

f_avg_unseg = transpose(mean(f_unseg, 2));

f_avg = transpose(mean(f, 2));
T_e = machine.R_s^2 * 2*pi * machine.L * imag(f_avg);

% Calculate zeroth radial order. Correct?
f_rr_unseg(1) = (max(real(f_avg_unseg)) - min(real(f_avg_unseg)))/2;
f_rr(1) = (max(real(f_avg)) - min(real(f_avg)))/2;
f_rt(1) = (max(imag(f_avg)) - min(imag(f_avg)))/2;

% Calculate cogging torque due to slots and segments
T_c_slot = cogging_torque(B_PM_sl(1,:), lambda_slot_total, machine, parts, theta_r, 'slot');
T_c_seg = cogging_torque(B_PM_sl(1,:), lambda_gap, machine, parts, theta_r, 'gap');
T_c = T_c_slot + T_c_seg;

% ------------------------------------------------------------------
%                        Mechanical calculations
% ------------------------------------------------------------------

% Calculate static deformations for modes r
Y_ms_unseg = deformation(machine, f_rr_unseg(1:25), r(1:25));
Y_ms = deformation(machine, f_rr(1:25), r(1:25));

f_n = Natural_frequencies(machine, r(1:25));

% ------------------------------------------------------------------
%                          Plot results
% ------------------------------------------------------------------

% Defaults
set(0,'defaulttextfontsize',20);
set(0,'defaultaxesfontsize',12);
set(0, 'DefaultLineLineWidth', 2);

set(0,'DefaultFigureColormap',parula);

theta_ticks = 0:2*pi/parts/4:2*pi/parts;
theta_labels = strrep(string(sym(theta_ticks)), 'pi', '\pi');

% Create colormap with white mapped to the lowest values
jet_copy = jet;
fft_colormap = [1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; jet_copy(20:end-20, :)];

% Show accuracy of Fourier series of permeance
[plt.lambda_re_slot, plt.lambda_im_slot] = ...
    Zarko_reconstruct(lambda_an_slot, lambda_bn_slot, machine.s, theta_vect_teeth, lambda_slot);
[plt.lambda_re_gap, plt.lambda_im_gap] = ...
    Zarko_reconstruct(lambda_an_gap, lambda_bn_gap, machine.Nseg, theta_vect, lambda_gap_c);

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
ylim([-0.4 0.4])
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_arm_sl_t = figure;
plot(theta_vect, imag(B_arm_sl(1,:)));
title('Slotless tangential flux density due to armature reaction')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
ylim([-0.4 0.4])
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_r_slotless = figure;
plot(theta_vect, real(B_slotless(1,:)));
title('Slotless radial flux density in full-load case')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot permeance 
plt.lambda = figure;
subplot(2,1,1);
plot(theta_vect, real(lambda));
title('Real part of the relative permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;
xlim([0 theta_ticks(end)])
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

subplot(2,1,2)
plot(theta_vect, imag(lambda));
title('Imaginary part of the relative permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;
xlim([0 theta_ticks(end)])
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)
plt.lambda.Position(3) = plt.lambda.Position(3)*2;

% Plot no-load flux density
plt.B_r_NL = figure;
plot(theta_vect, real(B_noload(1,:)));
title('Slotted radial flux density at t=0, no-load')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_t_NL = figure;
plot(theta_vect, imag(B_noload(1,:)));
title('Slotted tangential flux density at t=0, no-load')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

% Plot armature slotted flux density
plt.B_r_arm = figure;
plot(theta_vect, real(B_arm(1,:)));
title('Slotted radial flux density at t=0, armature')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_t_arm = figure;
plot(theta_vect, imag(B_arm(1,:)));
title('Slotted tangential flux density at t=0, armature')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_r_arm_fft = figure;
bar(k, [b_rk_arm; b_tk_arm], 1.5, 'FaceAlpha', 1);
set(gcf,'color','w');
title('Harmonics of the armature flux density')
xlabel('Spatial order')
ylabel('Magnetic flux density [T]')
legend('Radial', 'Tangential')
plt.B_r_arm_fft.Position(3) = plt.B_r_arm_fft.Position(3)*1.5;
ax=gca;
ax.XAxis.MinorTick = 'on';

% Plot on-load flux density
plt.B_r = figure;
plot(theta_vect, real(B(1,:)));
title('Slotted radial flux density at t=0')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_r_fft = figure;
bar(k, [b_rk_NL; b_rk], 1.5, 'FaceAlpha', 1);
set(gcf,'color','w');
title('Harmonics of the radial flux density')
xlabel('Spatial order')
ylabel('Magnetic flux density [T]')
legend('No-load', 'On-load')
plt.B_r_fft.Position(3) = plt.B_r_fft.Position(3)*1.5;
ax=gca;
ax.XAxis.MinorTick = 'on';

% Plot 2D FFT with decreased resolution for visibility
plt.B_r_fft_2d = figure;
imagesc(B_fx(1:2:end), B_ft(1:2:end), abs(B_2d_fft(1:2:end, 1:2:end)));
xlim([0, 150*parts])
ylim([-8 * machine.f_e, 8 * machine.f_e])
colormap(fft_colormap);
c = colorbar;
c.Label.String = 'Magnetic flux density [T]';
set(gcf,'color','w');
title('Harmonics of the radial flux density')
xlabel('Spatial order')
ylabel('Frequency [Hz]')

plt.B_t = figure;
plot(theta_vect, imag(B(1,:)));
title('Slotted tangential flux density at t=0')
xlabel('Angular position [rad]')
ylabel('Flux density [T]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.B_t_fft = figure;
bar(k, [b_tk_NL; b_tk], 1.5, 'FaceAlpha', 1);
set(gcf,'color','w');
title('Harmonics of the tangential flux density')
xlabel('Spatial order')
ylabel('Magnetic flux density [T]')
legend('No-load', 'On-load')
plt.B_t_fft.Position(3) = plt.B_t_fft.Position(3)*1.5;
ax=gca;
ax.XAxis.MinorTick = 'on';

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
set(gcf,'color','w');
title('Harmonics of the radial force')
xlabel('Spatial order')
ylabel('Force density [N/m^2]')
ax=gca;
ax.XAxis.MinorTick = 'on';

plt.f_r_fft_zoom = figure;
bar(r, [f_rr_unseg; f_rr], 'FaceAlpha', 1);
set(gcf,'color','w');
title('Harmonics of the radial force')
xlabel('Spatial order')
ylabel('Force density [N/m^2]')
xlim([-inf 320])
ylim([0 1500])
ax=gca;
ax.XAxis.MinorTick = 'on';
legend('Unsegmented', 'Segmented')
plt.f_r_fft_zoom.Position(3) = plt.f_r_fft_zoom.Position(3)*2;

plt.f_r_fft_2d = figure;
imagesc(f_fx(1:2:end), f_ft(1:2:end), abs(f_2d_fft(1:2:end, 1:2:end)));
xlim([0, 150*parts])
ylim([-8.2 * machine.f_e, 8.2 * machine.f_e])
colormap(fft_colormap);
c = colorbar;
c.Label.String = 'Force density [N/m^2]';
set(gcf,'color','w');
title('Harmonics of the radial force')
xlabel('Spatial order')
ylabel('Frequency [Hz]')

plt.f_t = figure;
plot(theta_vect, imag(f(1,:)));
title('Slotted tangential forces at t=0')
xlabel('Angular position [rad]')
ylabel('Force density [N/m^2]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

plt.f_t_fft = figure;
bar(r, f_rt, 'FaceAlpha', 1);
set(gcf,'color','w');
title('Harmonics of the tangential force')
xlabel('Spatial order')
ylabel('Force density [N/m^2]')
ax=gca;
ax.XAxis.MinorTick = 'on';

% Plot cogging torque for slots and slots+segments
plt.Tc = figure;
subplot(1,3,[1,2])
plot(theta_r, [T_c_slot; T_c]);
title('Cogging torque')
xlabel('Rotor angular position [rad]')
ylabel('Torque [Nm]')
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)

l_closeup = ceil(length(t_vect)/(machine.p*2/parts));   % Close up range
subplot(1,3,3)
plot(theta_r(1:l_closeup), [T_c_slot(1:l_closeup); T_c(1:l_closeup)]);
title('Close up')
xlabel('Rotor angular position [rad]')
ylabel('Torque [Nm]')
grid on;
legend({'Due to slots', "Due to slots"+newline+ "and segments"})
plt.Tc.Position(3) = plt.Tc.Position(3)*2;

% Plot average forces
plt.f_avg_r = figure;
plot(t_vect, real(f_avg));
title('Average radial force')
xlabel('Time [s]')
ylabel('Force density [N/m^2]')
ylim([0 max(real(f_avg))*1.1])
grid on;

axes('Position', [.58 .3 .3 .25])
box on
plot(t_vect(1:(length(f_avg)-1)/(machine.p/parts)), real(f_avg(1:(length(f_avg)-1)/(machine.p/parts))))
title(' One el. period')
xlabel('Time [s]')
ylabel('Force density [N/m^2]')

plt.Te_NL = figure;
plot(t_vect, T_e_NL);
title('Total torque, no-load case')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on;

plt.Te = figure;
plot(t_vect, T_e);
title('Total torque')
xlabel('Time [s]')
ylabel('Torque [Nm]')
ylim([min(T_e)*1.1 0])
grid on;

axes('Position', [.58 .5 .3 .25])
box on
plot(t_vect(1:(length(T_e)-1)/(machine.p/parts)), T_e(1:(length(T_e)-1)/(machine.p/parts)))
title(' One el. period')
xlabel('Time [s]')
ylabel('Torque [Nm]')

% Plot deformation
plt.Yms = figure;
bar(r(1:length(Y_ms)), [Y_ms_unseg; Y_ms], 'FaceAlpha', 1);
set(gcf,'color','w');
title('Quasi-static deformation for different modes')
xlabel('Spatial order')
ylabel('Deflection [m]')
legend('Unsegmented', 'Segmented')
plt.Yms.Position(3) = plt.Yms.Position(3)*1.5;

% Print natural frequencies
disp(f_n.Properties.Description)
disp(f_n)

% ------------------------------------------------------------------
%                            Save figures
% ------------------------------------------------------------------

% exportgraphics(plt.lambda_re_slot, 'Figures/lambda_re_slot.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.lambda_im_slot, 'Figures/lambda_im_slot.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.lambda_re_gap, 'Figures/lambda_re_gap.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.lambda_im_gap, 'Figures/lambda_im_gap.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_PM_sl, 'Figures/B_PM_sl.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_arm_sl_r, 'Figures/B_arm_sl_r.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_arm_sl_t, 'Figures/B_arm_sl_t.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_r_slotless, 'Figures/B_r_slotless.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.lambda, 'Figures/lambda.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_r, 'Figures/B_r.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_r_fft, 'Figures/B_r_fft.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_r_fft_2d, 'Figures/B_r_fft_2d.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_t, 'Figures/B_t.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.B_t_fft, 'Figures/B_t_fft.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_r, 'Figures/f_r.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_r_fft, 'Figures/f_r_fft.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_r_fft_zoom, 'Figures/f_r_fft_zoom.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_r_fft_2d, 'Figures/f_r_fft_2d.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_t, 'Figures/f_t.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_t_fft, 'Figures/f_t_fft.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.Tc, 'Figures/Tc.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.f_avg_r, 'Figures/f_avg_r.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.Te, 'Figures/Te.eps', 'BackgroundColor','none','ContentType','vector')
% exportgraphics(plt.Yms, 'Figures/Yms.eps', 'BackgroundColor','none','ContentType','vector')

% ------------------------------------------------------------------
%                             Create GIFs
% ------------------------------------------------------------------

%make_gif(theta_vect, t_vect, real(f), 'ForceRadial')