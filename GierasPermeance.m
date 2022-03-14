% Machine parameters
b14 = 0.002509;
g = 0.0005;
hM = 0.002;
t1 = 0.01;
muR = 1.05;

% Gieras method
kappa = b14/(g + hM/muR);
kappa2 = kappa/2;

gamma_1 = 4/pi * (kappa2 * atan(kappa2) - log(sqrt(1 + kappa2^2)));
kC = t1 / (t1 - gamma_1 * g);

gprime = g * kC + hM/muR;
rho = kappa/(kappa+5) * 2 * sqrt(1 + kappa^2) / (sqrt(1+kappa^2) - 1);

k = 1:20;
k_ok = sin(k*rho*pi*b14/(2*t1))./(k*rho*pi*b14/(2*t1));

Ak = -2 * gprime * gamma_1/t1 * k_ok.^2;

alpha_list = linspace(0, pi/18, 100);
lambda = ones(1, length(alpha_list));
for k_t = k
    lambda = lambda + Ak(k_t) * cos(k_t * 36 * (alpha_list + pi/36));
end

% Weber method
beta = 0.5 * (1 - 1/sqrt(1 + kappa^2));
Gamma = b14/t1;
Ak_Weber = -beta./k .* 4/pi .*(0.5 + (k*Gamma).^2 ./ (0.78 - 2*(k*Gamma).^2))...
    .* sin(1.6*pi*k*Gamma);

lambda_Weber = ones(1, length(alpha_list));
for k_t = k
    lambda_Weber = lambda_Weber + Ak_Weber(k_t) * cos(k_t * 36 * (alpha_list + pi/36));
end

figure
plot(alpha_list, lambda, alpha_list, lambda_Weber);