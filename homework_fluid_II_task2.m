%% Load Data 2
clear all;
close all;
atmosphere = load ('AtmosphericData_July_fs10Hz_Kurz.txt');
dataCenterline = load('Data_Centerline_FractalGrid_fs60kHz.txt');
load('fluctuations.mat');
%%
Fs_atmo = 10% sampling frequency atmo
Fs_dataCenter = 60% sampling frequency dataCenter

length_atmo = length(atmosphere);
length_dataCenter = length(dataCenterline);

fft_atmosphere = fft(atmosphere);
fft_dataCenterline = fft(dataCenterline);

p2_atmo = abs(fft_atmosphere);
p2_dataCenter = abs(fft_dataCenterline);

p1_atmo = p2_atmo(1:length_atmo/2+1);
p1_dataCenter = p2_dataCenter(1:length_dataCenter/2+1);

p1_atmo(2:end-1) = 2*p1_atmo(2:end-1);
p1_dataCenter(2:end-1) = 2*p1_dataCenter(2:end-1);

PSD_atmo = 2*p1_atmo.^2;
PSD_dataCenter = 2*p1_dataCenter.^2;

% testing for fluctuation
length_atmo_fluc = length(fluc_atmo);
length_center_fluc = length(fluc_center);

fft_atmosphere_fluc = fft(fluc_atmo);
fft_dataCenterline_fluc = fft(fluc_center);

p2_atmo_fluc = abs(fft_atmosphere_fluc);
p2_dataCenter_fluc = abs(fft_dataCenterline_fluc);

p1_atmo_fluc = p2_atmo_fluc(1:length_atmo_fluc/2+1);
p1_dataCenter_fluc = p2_dataCenter_fluc(1:length_center_fluc/2+1);

p1_atmo_fluc(2:end-1) = 2*p1_atmo_fluc(2:end-1);
p1_dataCenter_fluc(2:end-1) = 2*p1_dataCenter_fluc(2:end-1);

PSD_atmo_fluc = 2*p1_atmo_fluc.^2;
PSD_dataCenter_fluc = 2*p1_dataCenter_fluc.^2;


%%
f_atmo = Fs_atmo*(0:(length_atmo/2))/length_atmo;
figure;
plot(f_atmo,PSD_atmo)
title('One sided Amplitude Spectrum of u(t)')
xlabel('f (Hz)')
ylabel('|p1_atmo(f)|')

P_gerade_atmo = 10.^(-5/3 *log10(f_atmo)+4);
f_dataCenter = Fs_dataCenter*(0:(length_dataCenter/2))/length_dataCenter;
figure;
plot(f_dataCenter,PSD_dataCenter)
title('One sided Amplitude Spectrum of u(t)')
xlabel('f (Hz)')
ylabel('|p1_atmo(f)|')

P_gerade_dataCenter = 10.^(-5/3 *log10(f_dataCenter)+4);
figure
loglog(f_atmo,PSD_atmo)
hold on
plot(f_atmo,P_gerade_atmo)
title('PSD atmo')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo,2)])

figure
loglog(f_dataCenter,PSD_dataCenter)
hold on
plot(f_dataCenter,P_gerade_dataCenter)
title('PSD dataCenter')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_dataCenter,2)])

f_atmo_fluc = Fs_atmo*(0:(length_atmo_fluc/2))/length_atmo_fluc;
figure;
plot(f_atmo_fluc,PSD_atmo_fluc)
title('One sided Amplitude Spectrum of u_(t)')
xlabel('f (Hz)')
ylabel('|p1_atmo(f)|')

P_gerade_atmo_fluc = 10.^(-5/3 *log10(f_atmo_fluc)+4);
figure
loglog(f_atmo_fluc,PSD_atmo_fluc)
hold on
plot(f_atmo_fluc,P_gerade_atmo_fluc)
title('PSD atmo fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo_fluc,2)])

f_center_fluc = Fs_dataCenter*(0:(length_center_fluc/2))/length_center_fluc;
figure;
plot(f_center_fluc,PSD_dataCenter_fluc)
title('One sided Amplitude Spectrum of u_(t)')
xlabel('f (Hz)')
ylabel('|p1_center(f)|')

P_gerade_center_fluc = 10.^(-5/3 *log10(f_center_fluc)+4);
figure
loglog(f_center_fluc,PSD_dataCenter_fluc)
hold on
plot(f_center_fluc,P_gerade_center_fluc)
title('PSD Center fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center_fluc,2)])

%% Smoothing
% measured velocity atmosphere
PSD_atmo_smooth = smooth(PSD_atmo,20);
figure
loglog(f_atmo,PSD_atmo_smooth)
hold on
plot(f_atmo,P_gerade_atmo)
title('smoothed PSD atmosphere')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo,2)])
% fluctuation atmosphere

PSD_atmo_fluc_smooth = smooth(PSD_atmo_fluc,20);
figure
loglog(f_atmo_fluc,PSD_atmo_fluc_smooth)
hold on
plot(f_atmo_fluc,P_gerade_atmo_fluc)
title('smoothed PSD atmo fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo_fluc,2)])

%% Statistics
[x,y] = ksdensity(fluc_atmo);
%% Correlation
%integral length
timeLags = 12500;
autocorr_time_lag = (1-1:timeLags)';
autocorr_data = autocorr(fluc_atmo,timeLags);

% check fft of autocorr
autocorr_fft = fft(autocorr_data);
figure
plot(abs(autocorr_fft(1:length(autocorr_fft)/2+1)));
%delete data lower than 0 of the autocorrelation (visually dicided)
i = 1
while autocorr_data(i,1) > 0
    i = i + 1
end
autocorr_time_lag(i:end) = [];
autocorr_data(i:end) = [];
figure
plot(autocorr_data);
%% Joint probability dist
xlag = lagmatrix(fluc_atmo,[0 1 2 3 4 5 6 7 8 9 10]);

x_axis = -3:.2:3; % Define edges of bins for x axis. Column vector
y_axis = -3:.2:3; % Same for y axis

%// Compute and plot pdf
figure
histogram2(xlag(:,1), xlag(:,2), x_axis, y_axis, 'Normalization', 'pdf')

%// Compute and plot pdf
figure
histogram2(xlag(:,1), xlag(:,10), x_axis, y_axis, 'Normalization', 'pdf')

%// Compute and plot cdf
figure
histogram2(xlag(:,1), xlag(:,2), x_axis, y_axis, 'Normalization', 'cdf')
%%
% hist3 will bin the data
xi = linspace(min(xlag(:,1)), max(xlag(:,2)), 50);
yi = linspace(min(xlag(:,1)), max(xlag(:,2)), 50);
hst = hist3(xlag(:,1:2),{xi yi}); %removed extra '

% normalize the histogram data
dx = xi(2)-xi(1);
dy = yi(2)-yi(1);
area = dx*dy;
pdfData = hst/sum(sum(hst))/area;

% plot pdf
figure; 
contour(xi,yi,pdfData);

vector(:,1) = xlag(:,1);
vector(:,2) = xlag(:,10)
% hist3 will bin the data
xi = linspace(min(xlag(:,1)), max(xlag(:,10)), 50);
yi = linspace(min(xlag(:,1)), max(xlag(:,10)), 50);
hst = hist3(vector(:,1:2),{xi yi}); %removed extra '

% normalize the histogram data
dx = xi(2)-xi(1);
dy = yi(2)-yi(1);
area = dx*dy;
pdfData = hst/sum(sum(hst))/area;

% plot pdf
figure; 
contour(xi,yi,pdfData);

