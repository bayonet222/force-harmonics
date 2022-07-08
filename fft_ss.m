function [m, F_m] = fft_ss(signal, N_F, multiplier)
% Single-sided FFT of the signal with N_F Fourier components for ease of using fft
%   Multiplier: amount of times signal repeats

L = length(signal);

% Define single-sided FFT
Y = fft(signal);
P2 = abs(Y/ L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

% Output FFT and orders until N_F
F_m = P1(1:N_F);
m = (0:N_F-1)*multiplier;
end

