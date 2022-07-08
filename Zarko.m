function [lambda, lambda_an, lambda_bn] = Zarko(machine_params, theta_list, N_lambda, type)
% Function to analyse the relative permeance due to slots and gaps

% Load machine paramaters
Rs = machine_params.R_s;                    % Stator inner radius
Rr = machine_params.R_r;                    % Rotor outer radius

% Analyse either slots or segment gaps
if isequal(type, 'slot')
    Ns = machine_params.s;                  % Number of slots
    theta_o = machine_params.theta_o;       % Slot opening pitch
elseif isequal(type, 'gap')
    Ns = machine_params.Nseg;               % Number of segments
    theta_o = machine_params.theta_gap;     % Gap width angle
end

% ---------- Zarko -----------
theta_s = 2*pi/Ns;                          % Pitch
theta_1 = theta_s/2 - theta_o/2;            % Inner slot end
theta_2 = theta_s/2 + theta_o/2;            % Outer slot end

g = log(Rs/Rr);
R_eval = Rs - (Rs - Rr)/50;         % Evaluate lambda just below surface
b_o = theta_o;                      % Opening pitch
b = (b_o/(2*g) + sqrt((b_o/(2*g))^2 + 1))^2;
a = 1/b;

C = log(Rs) + 1j * theta_2;

% Define equation to be solved numerically
syms w
p = sqrt((w-b)/(w-a));
z = 1j * g/pi * (log((1+p)/(1-p)) - log((b+p)/(b-p)) - 2*(b-1)/sqrt(b) * atan(p/sqrt(b))) + C;

% Define 'interesting region' to speed up calculations
if isequal(type, 'slot')
    region_b = 0;
    region_e = theta_s;
elseif isequal(type, 'gap')
    region_b = theta_s/2 - 3 * theta_o;
    region_e = theta_s/2 + 3 * theta_o;
end

% Most of lambda equals 1 + 0j
lambda = ones(1, length(theta_list)) * (1 + 0j);
samples = length(theta_list);

% Solve Zarko's formula for each value of theta in the region of interest
for n = 1:length(theta_list)
    theta = theta_list(n);
    if (region_b <= theta) && (theta <= region_e)
        s = R_eval * exp(1j * theta);
        w_solved = vpasolve(z == log(s), w);
        w_solved = double(w_solved);

        if isempty(w_solved)
            lambda(n) = 1.0 + 0.0j;     % Assume lambda is 1 if unsolvable
        else
            k = Rs * exp(1j * (g/pi * log(w_solved) + theta_s/2));
            lambda(n) = k/s * (w_solved-1)/(sqrt(w_solved-a)*sqrt(w_solved-b));
        end
    end
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

lambda_an = real(P1_re(1:N_lambda+1));
lambda_bn = imag(P1_im(1:N_lambda+1));
end