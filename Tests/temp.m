% Fourier parameters
samples = 101;
N_lambda = floor(samples/2); % Number of Fourier terms to consider, max
N_lambda = 20;

% Geometric parameters
Rs = 0.0575; 
Rr = 0.055;

Qs = 36;
theta_s = 2*pi/Qs;
theta_o = theta_s/4;    % Angular slot opening

% ---------- Zarko -----------
theta_1 = theta_s/2 - theta_o/2;
theta_2 = theta_s/2 + theta_o/2;

g = log(Rs/Rr);
b_o = theta_o;
b = (b_o/(2*g) + sqrt((b_o/(2*g))^2 + 1))^2;
a = 1/b;

C = log(Rs) + 1j * theta_2;

syms w
p = sqrt((w-b)/(w-a));
z = 1j * g/pi * (log((1+p)/(1-p)) - log((b+p)/(b-p)) - 2*(b-1)/sqrt(b) * atan(p/sqrt(b))) + C;

lambda = [];
theta_list = linspace(0, theta_s, samples);

% Solve Zarko's formula for each value of theta
for theta = theta_list
    s = (Rs - 0.0001) * exp(1j * theta);
    w_solved = eval(vpasolve(z == log(s)));

    k = Rs * exp(1j * (g/pi * log(w_solved) + theta_s/2));
    lambda_new = k/s * (w_solved-1)/(sqrt(w_solved-a)*sqrt(w_solved-b));
    lambda = [lambda, lambda_new];
end

% Fourier transform
FFT_re = fft(real(lambda.'));
FFT_im = fft(imag(lambda.'));

% Real part
P2_re = FFT_re/samples;
P1_re = P2_re(1:floor(samples/2)+1);
P1_re(2:end) = 2 * P1_re(2:end);

% Imaginary part
P2_im = FFT_im/samples;
P1_im = P2_im(1:floor(samples/2)+1);
P1_im(2:end) = 2 * P1_im(2:end);

% Plot frequency spectrum
f = 1/theta_s * (0:floor(samples/2));
figure
plot(f, real(P1_re))
title('Frequency spectrum of real part')
figure
plot(f, imag(P1_im))
title('Frequency spectrum of imaginary part')

% Initialize reconstruction of lambda with a0
% lambda_rec_re = ones(1, samples) * real(P2_re(1));
% lambda_rec_im = zeros(1,samples);

% Reconstruct lambda to verify if the calculation was correct
% for i = 1:N_lambda
%     lambda_rec_re = lambda_rec_re + real(P1_re(i+1)) * cos(i * Qs * theta_list) + 1j * imag(P1_re(i+1)) * sin(i * Qs * theta_list);
%     lambda_rec_im = lambda_rec_im + real(P1_im(i+1)) * cos(i * Qs * theta_list) + 1j * imag(P1_im(i+1)) * sin(i * Qs * theta_list);
% end

% Take the real component of the real fft and the negative imaginary
% component of the imag fft. Not 100% sure why but it works
% lambda_rec = real(lambda_rec_re) - 1j * imag(lambda_rec_im);

lambda_an = real(P1_re(1:N_lambda+1));
lambda_bn = imag(P1_im(1:N_lambda+1));

% Could be simplified to
lambda_rec = ones(1,samples) * lambda_an(1);
for i = 1:N_lambda
    lambda_rec = lambda_rec + lambda_an(i+1) * cos(i * Qs * theta_list) - 1j * lambda_bn(i+1) * sin(i * Qs * theta_list);
end

% Plot the reconstructed lambda
figure
plot(theta_list, real(lambda), theta_list, real(lambda_rec))
legend('Exact', 'Reconstructed')
title('Real part of the complex permeance')

figure
plot(theta_list, imag(lambda), theta_list, imag(lambda_rec))
legend('Exact', 'Reconstructed')
title('Imaginary part of the complex permeance')