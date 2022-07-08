function [handle_re, handle_im] = Zarko_reconstruct(lambda_an, lambda_bn, Ns, theta_list, lambda)
% Create plot showing FFT and original function for lambda

% Multiply i by Ns as the waveform repeats Ns times along the
% circumference of the machine
lambda_rec = ones(1,length(theta_list)) * lambda_an(1);
for i = 1:length(lambda_an)-1
    lambda_rec = lambda_rec + lambda_an(i+1) * cos(i * Ns * theta_list) - 1j * lambda_bn(i+1) * sin(i * Ns * theta_list);
end

% Plot the reconstructed lambda
handle_re = figure;
plot(theta_list, real(lambda), theta_list, real(lambda_rec), '--')
legend('Exact', 'Fourier series')
title('Real part of the complex permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;

handle_im = figure;
plot(theta_list, imag(lambda), theta_list, imag(lambda_rec), '--')
legend('Exact', 'Fourier series')
title('Imaginary part of the complex permeance')
xlabel('Angular position [rad]')
ylabel('Relative permeance')
grid on;
end