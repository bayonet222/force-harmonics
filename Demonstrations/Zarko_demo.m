% THIS FILE SHOWS HOW TO USE THE ZARKO AND ZARKO-RECONSTRUCT FILES

% Fourier parameters
samples = 101;
N_lambda = floor(samples/2); % Number of Fourier terms to consider, max
N_lambda = 20;

% Geometric parameters
machine_params.Rs = 0.0575; 
machine_params.Rr = 0.055;

machine_params.s = 36;
theta_s = 2*pi/machine_params.s;
theta_list = linspace(0, theta_s, samples);
machine_params.theta_o = theta_s/4;    % Angular slot opening

[lambda, lambda_an, lambda_bn] = Zarko(machine_params, theta_list, N_lambda);

Zarko_reconstruct(lambda_an, lambda_bn, machine_params.s, theta_list, lambda)