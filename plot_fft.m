function plot_fft(signal, N_F, multiplier)
%Barplot of the Fast Fourier Transform of the signal

[m, F_m] = fft_ss(signal, N_F, multiplier);
bar(m, F_m)
end

