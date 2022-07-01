function export_force(theta_vect, t_vect, f)
% Define time vector length and period
L = length(t_vect);
T = t_vect(end);

% Frequency range
freq_vect = [0 : floor(L/2)]/ T;

% FFT of radial force
Yr = fft(real(f));
P2r = Yr/ L;
P1r = P2r(1:floor(L/2)+1, :);
P1r(2:end-1, :) = 2 * P1r(2:end-1, :);

% Create 1D vector, stacking columns of theta
Pr_list = P1r(:);

% FFT of radial force
Yt = fft(imag(f));
P2t = Yt/ L;
P1t = P2t(1:floor(L/2)+1, :);
P1t(2:end-1, :) = 2 * P1t(2:end-1, :);

% Create 1D vector, stacking columns of theta
Pt_list = P1t(:);

% Create combinations of theta and frequencies
[a, b] = meshgrid(theta_vect, freq_vect);

% Export theta, frequency, force amplitude and force phase (radial and
% tangential)
out = [a(:) b(:) abs(Pr_list) angle(Pr_list) abs(Pt_list) angle(Pt_list)];

writematrix(out, 'forces.txt')
writematrix(freq_vect, 'freqs.txt')

end