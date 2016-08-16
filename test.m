x = load ('AtmosphericData_July_fs10Hz_Kurz.txt');

Fs = 10

n = length(atmosphere)

% cross correlation reference
xref = xcorr(x,x);

%FFT method based on zero padding
fx = fft([x; zeros(n,1)]); % zero pad and FFT
x2 = ifft(fx.*conj(fx)); % abs()^2 and IFFT
% circulate to get the peak in the middle and drop one
% excess zero to get to 2*n-1 samples
x3 = [x2(n+2:end); x2(1:n)];
% calculate the error
d = x3-xref; % difference, this is actually zero
fprintf('Max error = %6.2f\n',max(abs(d)));

figure
plot(x2)
figure
plot(x)
