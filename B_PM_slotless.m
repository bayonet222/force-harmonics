
function [B_mmu, B_PM] = B_PM_slotless(theta_vect, t_vect, machine_params)
%Function B_PM of theta
%   Returns the function of the magnetic flux density
%   due to the permanent magnets at the stator teeth
%   as function of theta


% Load machine parameters
p = machine_params.p;
R_s = machine_params.R_s;
R_r = machine_params.R_r;
R_m = machine_params.R_m;
Br = machine_params.Br;
alpha_m = machine_params.alpha_m;
mu_r = machine_params.mu_r;

omega_m = machine_params.omega_m;       % Rotation speed

% Set mu
mu_max = 35;
mu = transpose(1:2:mu_max);            % Odd integers

B_mmu = K_B(mu, p, R_s, R_r, R_m, Br, alpha_m, mu_r) ...
    .* f_Br(R_s, R_m, mu, p); 

% Obsolete
B_PM_t0 = sum(B_mmu .* cos(mu .* p .* theta_vect));

% Create theta-t matrix and loop through time
B_PM = zeros(length(t_vect), length(theta_vect));

for i = 1:length(t_vect)
    B_PM(i,:) = sum(B_mmu .* cos(mu.*p.*(theta_vect - omega_m*t_vect(i))));
end

end

function K_B_mu = K_B(mu, p, R_s, R_r, R_m, Br, alpha_m, mu_r)
global mu_0

K_B_mu = mu_0/mu_r * M_mu(mu, Br, alpha_m) .* mu*p./((mu*p).^2 - 1) .* ...
    (mu*p - 1 + 2*(R_r/R_m).^(mu*p+1) - (mu*p+1).*(R_r/R_m).^(2*mu*p))./...
    ((mu_r +1)/mu_r * (1-(R_r/R_s).^(2*mu*p)) - ...
        (mu_r-1)/mu_r * ((R_m/R_s).^(2*mu*p) - (R_r/R_m).^(2*mu*p)));
end

function magnetization = M_mu(mu, Br, alpha_m)
global mu_0

magnetization = 4*Br/mu_0 * alpha_m *sin(mu * pi * alpha_m/2)./...
    (mu * pi * alpha_m);

end

function f_Br_Rs = f_Br(R_s, R_m, mu, p)
f_Br_Rs = 2 * (R_m/R_s).^(mu*p + 1);
end