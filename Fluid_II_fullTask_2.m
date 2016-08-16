%% two-point quantities
% load data
clear all;
close all;
atmosphere = load ('AtmosphericData_July_fs10Hz_Kurz.txt');
dataCenterline = load('Data_Centerline_FractalGrid_fs60kHz.txt');
load('fluctuations.mat');
%% power spectrum E(f) or E(k)
Fs_atmo = 10    % sampling frequency atmo

%   atmosphere u
length_atmo = length(atmosphere);
fft_atmosphere = fft(atmosphere);
p2_atmo = abs(fft_atmosphere);
p1_atmo = p2_atmo(1:length_atmo/2+1);
p1_atmo(2:end-1) = 2*p1_atmo(2:end-1);
PSD_atmo = 2*p1_atmo.^2;

f_atmo = Fs_atmo*(0:(length_atmo/2))/length_atmo;
figure;
plot(f_atmo,PSD_atmo)
title('Atmosphere - One sided Amplitude Spectrum of u(t)')
xlabel('f (Hz)')
ylabel('|p1_atmo(f)|')

%   atmospere fluc u
length_atmo_fluc = length(fluc_atmo);
fft_atmosphere_fluc = fft(fluc_atmo);
p2_atmo_fluc = abs(fft_atmosphere_fluc);
p1_atmo_fluc = p2_atmo_fluc(1:length_atmo_fluc/2+1);
p1_atmo_fluc(2:end-1) = 2*p1_atmo_fluc(2:end-1);
PSD_atmo_fluc = 2*p1_atmo_fluc.^2;

f_atmo_fluc = Fs_atmo*(0:(length_atmo_fluc/2))/length_atmo_fluc;
figure;
plot(f_atmo_fluc,PSD_atmo_fluc)
title('Atmosphere - One sided Amplitude Spectrum of u dash (t)')
xlabel('f (Hz)')
ylabel('|p1_atmo(f)|')

Fs_center = 60% sampling frequency dataCenter

%   dataCenter u
length_center = length(dataCenterline);
fft_center = fft(dataCenterline);
p2_center = abs(fft_center);
p1_center = p2_center(1:length_center/2+1);
p1_center(2:end-1) = 2*p1_center(2:end-1);
PSD_center = 2*p1_center.^2;

f_center = Fs_center*(0:(length_center/2))/length_center;
figure;
plot(f_center,PSD_center)
title('dataCenter - One sided Amplitude Spectrum of u(t)')
xlabel('f (Hz)')
ylabel('|p1_center(f)|')

%   dataCenter fluc u
length_center_fluc = length(fluc_center);
fft_center_fluc = fft(fluc_center);
p2_center_fluc = abs(fft_center_fluc);
p1_center_fluc = p2_center_fluc(1:length_center_fluc/2+1);
p1_center_fluc(2:end-1) = 2*p1_center_fluc(2:end-1);
PSD_center_fluc = 2*p1_center_fluc.^2;

f_center_fluc = Fs_center*(0:(length_center_fluc/2))/length_center_fluc;
figure;
plot(f_center_fluc,PSD_center_fluc)
title('dataCenter - One sided Amplitude Spectrum of u dash (t)')
xlabel('f (Hz)')
ylabel('|p1_center(f)|')

%% fit for 5/3
P_gerade_atmo = 10.^(-5/3 *log10(f_atmo)+4)
P_gerade_atmo_fluc = 10.^(-5/3 *log10(f_atmo_fluc)+4);
P_gerade_center = 10.^(-5/3 *log10(f_center)+4);
P_gerade_center_fluc = 10.^(-5/3 *log10(f_center_fluc)+4);

% fit plot atmosphere u
figure
loglog(f_atmo',PSD_atmo)
hold on
plot(f_atmo',P_gerade_atmo)
title('PSD atmo')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo',1)])

% fit plot atmosphere fluc u 
figure
loglog(f_atmo_fluc',PSD_atmo_fluc)
hold on
plot(f_atmo_fluc',P_gerade_atmo_fluc)
title('PSD atmo fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo_fluc',1)])

% fit for plot dataCenter u
figure
loglog(f_center',PSD_center)
hold on
plot(f_center',P_gerade_center)
title('PSD dataCenter')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center',1)])

% fit for plot dataCenter fluc u
figure
loglog(f_center_fluc',PSD_center_fluc)
hold on
plot(f_center_fluc',P_gerade_center_fluc)
title('PSD dataCenter fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center_fluc',1)])

%% fft of autocorrelation

%atmosphere
autocorr_atmo  = xcorr(atmosphere,atmosphere);
fft_auto_corr = fft(autocorr_atmo);

figure
loglog(abs(fft_auto_corr));
disp('test')
hold on
loglog(PSD_atmo) % DONE here

Y = fft(atmosphere);
L = length(atmosphere);
P2 = abs(Y);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
PSD = 2*P1.^2;

rev_PSD = ifft(PSD);
figure
plot(1:length(rev_PSD),abs(rev_PSD));
