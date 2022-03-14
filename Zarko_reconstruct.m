function [handle_re, handle_im] = Zarko_reconstruct(lambda_an, lambda_bn, Ns, theta_list, lambda)
    % Old method to reconstruct full real and imaginary components
    % Initialize reconstruction of lambda with a0
    % lambda_rec_re = ones(1, samples) * real(P2_re(1));
    % lambda_rec_im = zeros(1,samples);

    % Reconstruct lambda to verify if the calculation was correct
    % for i = 1:N_lambda
    %     lambda_rec_re = lambda_rec_re + real(P1_re(i+1)) * cos(i * Ns * theta_list) + 1j * imag(P1_re(i+1)) * sin(i * Ns * theta_list);
    %     lambda_rec_im = lambda_rec_im + real(P1_im(i+1)) * cos(i * Ns * theta_list) + 1j * imag(P1_im(i+1)) * sin(i * Ns * theta_list);
    % end

    % Take the real component of the real fft and the negative imaginary
    % component of the imag fft. Not 100% sure why but it works
    % lambda_rec = real(lambda_rec_re) - 1j * imag(lambda_rec_im);
    
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