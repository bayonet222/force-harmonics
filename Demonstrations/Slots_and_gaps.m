% THIS FILE SHOWS HOW TO COMBINE SLOTS AND FLUX GAPS

% Fourier parameters
samples = 200;
N_lambda = floor(samples/2); % Number of Fourier terms to consider, max
N_lambda = 100;

% Geometric parameters
Rs = 0.0575; 
Rr = 0.055;

Ns = 36;
theta_s = 2*pi/Ns;
theta_o = theta_s/4;    % Angular slot opening

Ng = 9;
theta_seg = 2 * pi/Ng;
theta_og = 0.01;        % Angular flux gap opening

[lambda_slot, theta_list_slot, lambda_an_slot, lambda_bn_slot] = Zarko(Rs, Rr, Ns, theta_o, samples, N_lambda);

% Zarko_reconstruct(lambda_an_slot, lambda_bn_slot, Ns, theta_list_slot, lambda_slot)

[lambda_gap, theta_list_gap, lambda_an_gap, lambda_bn_gap] = Zarko(Rs, Rr, Ng, theta_og, samples, N_lambda);

% Zarko_reconstruct(lambda_an_gap, lambda_bn_gap, Ng, theta_list_gap, lambda_gap)

lambda_rec = ones(1,length(theta_list_gap)) * lambda_an_slot(1);
for i = 1:N_lambda
    lambda_rec = lambda_rec + lambda_an_slot(i+1) * cos(i * Ns * theta_list_gap) ...
        - 1j * lambda_bn_slot(i+1) * sin(i * Ns * theta_list_gap)...
        + lambda_an_gap(i+1) * cos(i * Ng * theta_list_gap)...
        - 1j * lambda_bn_gap(i+1) * sin(i * Ng * theta_list_gap);
end

% Interpolate the repeated slot permeance to match the theta list
% The gaps repeat Ns/Ng times over 1 segment, and overlap at the begin and
% end point of each slot span
lambda_slot_interp = interp1(linspace(0, Ns/Ng * theta_s, Ns/Ng * (samples - 1) + 1), [repmat(lambda_slot(1:end-1), 1, Ns/Ng), lambda_slot(end)], theta_list_gap);

% Plot the reconstructed lambda
figure
plot(theta_list_gap, real(lambda_slot_interp + lambda_gap) - 1, theta_list_gap, real(lambda_rec))
legend('Exact', 'Reconstructed')
title('Real part of the complex permeance')

figure
plot(theta_list_gap, imag(lambda_slot_interp + lambda_gap), theta_list_gap, imag(lambda_rec))
legend('Exact', 'Reconstructed')
title('Imaginary part of the complex permeance')