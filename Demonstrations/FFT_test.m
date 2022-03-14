clear;
clc;
Fs = 1000;          %sampling rate
T = 1/Fs;           %sampling interval
L = 500;            %Number of time points
t = 0:T:(L-1)*T;    %time vector
N = 1024;           %FFT points
% Define frequencies
f1 = 20;
f2 = 220;
f3 = 138; 
% Create the signal
x  = 5 + 12*sin(2*pi*f1*t-0.3) + 5*sin(2*pi*f2*t-0.5) + 4*sin(2*pi*f3*t - 1.5);
figure(1)
plot(t*1000,x)
xlabel('time (in ms)')
ylabel('Amplitude')
title('Original Signal')
% Add some noise
figure(2)
% x = awgn(x,5,'measured')
x = x + 2 * randn(size(t));
plot(t*1000,x)
xlabel('time (in ms)')
ylabel('Amplitude')
title('Noisy Signal')

% FFT
X = fft(x,N);
figure(3)
plot(abs(X))
SSB = X(1:N/2);
SSB(2:end) = 2*SSB(2:end);
f = (0:N/2-1)*(Fs/N);
% Amplitude
figure(4)
plot(f,abs(SSB/L))
xlabel('f (in Hz)')
ylabel('|X(f)|')
title('Corrected Frequency Spectrum')