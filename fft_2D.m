function [fx, ft, Y2] = fft_2D(Data,Xmultiplier, T)
%FFT_2D Summary of this function goes here
%   X multiplier: number of times X repeats
%   T: period of the time signal

[lenT, lenX] = size(Data);
Y = fftshift(fft2(Data))/numel(Data);
% Create single-sided spectrum
Y2 = Y + Y(end:-1:1, end:-1:1);

% Remove DC component
Y2(ceil(lenT/2), ceil(lenX/2)) = 0;

fx = ([0 : lenX-1] - (lenX - 1)/2) * Xmultiplier;
ft = ([0 : lenT-1] - (lenT - 1)/2) / T;

end

