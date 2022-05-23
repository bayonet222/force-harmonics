function T_c = cogging_torque(B_PM, lambda, machine_params, parts, theta_r, type)
% Function to calculate the cogging torque

    global mu_0

    % Machine params
    p = machine_params.p;
    L = machine_params.L;
    g = machine_params.g;
    R_m = machine_params.R_m;
    
    if isequal(type, 'slot')
        s = machine_params.s;
        theta_c = 0;
    elseif isequal(type, 'gap')
        s = machine_params.Nseg;
        s2 = machine_params.s;
        theta_c = -pi/lcm(s2, 2*p)/4;       % phase shift for gap CT
    end
    
    N_L = lcm(s, 2*p);

    % Obtain the FFT of the squared flux and permeance
    [~, B_mmu2] = fft_ss(B_PM.^2, floor(length(B_PM)/2), 1);
    [~, lambda_an2] = fft_ss(real(lambda).^2, floor(length(lambda)/2), 1);
    [~, lambda_n2] = fft_ss(real(lambda).*imag(lambda), floor(length(lambda)/2), 1);

    B_nNL = transpose(B_mmu2(N_L/parts + 1:N_L/parts:end));
    lambda_nNL = transpose(lambda_an2(N_L/parts + 1:N_L/parts:end));
    lambda_tnNL = transpose(lambda_n2(N_L/parts + 1:N_L/parts:end));

    N_CT = min(length(B_nNL), length(lambda_nNL));  % Number of CT harmonics
    kNL = transpose(1:N_CT)*N_L;                    % Cogging torque harmonics

    T_c = L * R_m * g * pi / 2 / mu_0 * ...
        sum(kNL .* B_nNL(1:N_CT) .* lambda_nNL(1:N_CT) .* sin(kNL .* (theta_r + theta_c)));
    
    T_c = L * R_m^2 * pi/ mu_0 * ...
        sum(B_nNL(1:N_CT) .* lambda_tnNL(1:N_CT) .* sin(kNL .* (theta_r + theta_c)));
end