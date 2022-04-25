function [fx, ft, Y] = fft_2D(Data,Xmultiplier, T)
%FFT_2D Summary of this function goes here
%   X multiplier: number of times X repeats
%   T: period of the time signal

[lenT, lenX] = size(Data);
Y = fftshift(fft2(Data))/numel(Data);

fx = ([0 : lenX-1] - (lenX - 1)/2) * Xmultiplier;
ft = ([0 : lenT-1] - (lenT - 1)/2) / T;

end

