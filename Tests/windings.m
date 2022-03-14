s = 12;
beta_s = 2*pi/s;

theta = linspace(0, 2*pi, s*100);
n_phi = zeros(1, length(theta));

for t = [1 650 1150]
    n_phi(t) = 1;
end

for t = [50 550]
    n_phi(t) = -1;
end

N_c = cumtrapz(theta, n_phi)*s*100/2/pi;

N = N_c - (sum(N_c)/length(N_c));

% FFT
L = length(N);

Y = fft(N);
P2 = Y/ L;
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

N_F = 50;

k_w = P1(2:N_F+1);
% This is equal to sin(v * beta/2)/v * 4/pi

% Get FFT of full pitch
h = 1:N_F;
k_w_FP = zeros(1, N_F);

for i = h
    if mod(i, 2) == 1
        k_w_FP(i) = (-1)^((i-1)/2) * 4/i/pi;
    end
end
% This is equal to +/- 1/v * 4/pi

figure
bar(theta, n_phi)

figure
plot(theta, N_c)

figure
plot(theta, N)

figure
bar(h, k_w./k_w_FP)