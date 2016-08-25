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
P_gerade_atmo = 10.^(-5/3 *log10(f_atmo)+4);
P_gerade_atmo_fluc = 10.^(-5/3 *log10(f_atmo_fluc)+4);
P_gerade_center = 10.^(-5/3 *log10(f_center)+4);
P_gerade_center_fluc = 10.^(-5/3 *log10(f_center_fluc)+4);

% fit plot atmosphere u
figure
loglog(f_atmo',PSD_atmo);
hold on
plot(f_atmo',P_gerade_atmo)
title('PSD atmo')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo',1)])

% fit plot atmosphere fluc u 
figure
loglog(f_atmo_fluc',PSD_atmo_fluc);
hold on
plot(f_atmo_fluc',P_gerade_atmo_fluc)
title('PSD atmo fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo_fluc',1)])

% fit for plot dataCenter u
figure
loglog(f_center',PSD_center);
hold on
plot(f_center',P_gerade_center)
title('PSD dataCenter')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center',1)])

% fit for plot dataCenter fluc u
figure
loglog(f_center_fluc',PSD_center_fluc);
hold on
plot(f_center_fluc',P_gerade_center_fluc)
title('PSD dataCenter fluc')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center_fluc',1)])

%% fft of autocorrelation, numerically equal

%atmosphere
autocorr_atmo  = xcorr(atmosphere,atmosphere,length(f_atmo));
fft_auto_corr_atmo = fft(autocorr_atmo);
fft_auto_corr_atmo = abs(fft_auto_corr_atmo);
fft_auto_corr_atmo = fft_auto_corr_atmo(1:length(atmosphere)/2+1);

figure
loglog(f_atmo,fft_auto_corr_atmo);
hold on
loglog(f_atmo,PSD_atmo); % DONE here
title('PSD of atmosphere')
xlabel('f (Hz)')
ylabel('Power spectral density')
legend('FFT of autocorrelation of atmosphere','FFT squared');

%dataCenter
autocorr_center  = xcorr(dataCenterline,dataCenterline);
fft_auto_corr_center = fft(autocorr_center);
fft_auto_corr_center = abs(fft_auto_corr_center);
fft_auto_corr_center = fft_auto_corr_center(1:length(dataCenterline)/2+1);

figure
loglog(f_center,fft_auto_corr_center);
hold on
loglog(f_center,PSD_center); % DONE here
title('PSD of dataCenterline')
xlabel('f (Hz)')
ylabel('Power spectral density')
legend('FFT of autocorrelation of dataCenterline','FFT squared');
%%  Joint probability dist
% atmosphere , create lagmatrix
xlag_atmo = lagmatrix(fluc_atmo(:,1),[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);
% check for correlation

figure
for i = 1 : length(xlag_atmo(1,:))
    vector(:,1) = xlag_atmo(:,1);
    vector(:,2) = xlag_atmo(:,i);
    xi = linspace(min(vector(:,1)), max(vector(:,1)), 50);
    yi = linspace(min(vector(:,2)), max(vector(:,2)), 50);
    hst = hist3(vector(:,1:2),{xi yi}); %removed extra '
    %normalize
    dx = xi(2)-xi(1);
    dy = yi(2)-yi(1);
    area = dx*dy;
    pdfData = hst/sum(sum(hst))/area;
    subplot(4,4,i)
    contour(xi,yi,pdfData);
    title(['shift of' num2str(i)])
end

% dataCenterline , create lagmatrix
xlag_dataCenter = lagmatrix(fluc_center(:,1),[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]);
% check for correlation

figure
for i = 1 : length(xlag_dataCenter(1,:))
    vector_center(:,1) = xlag_dataCenter(:,1);
    vector_center(:,2) = xlag_dataCenter(:,i);
    xi = linspace(min(vector_center(:,1)), max(vector_center(:,1)), 50);
    yi = linspace(min(vector_center(:,2)), max(vector_center(:,2)), 50);
    hst = hist3(vector_center(:,1:2),{xi yi}); %removed extra '
    %normalize
    dx = xi(2)-xi(1);
    dy = yi(2)-yi(1);
    area = dx*dy;
    pdfData = hst/sum(sum(hst))/area;
    subplot(4,4,i)
    contour(xi,yi,pdfData);
    title(['shift of' num2str(i)])
end


%% Taylor's hypothesis

r_atmo = nanmean(atmosphere) / Fs_atmo; % diameter of frozen turbulence structure
r_dataCenter = nanmean(dataCenterline) / Fs_center; % diameter of frozen turbulence structure
%% Integral Length
%   atmo
timeLags_atmo = length(atmosphere)-1;
autocorr_time_lag_atmo = (1-1:timeLags_atmo)';
autocorr_data_atmo = autocorr(atmosphere,timeLags_atmo);
%delete data lower than 0 of the autocorrelation (visually dicided)
i = 1;
while autocorr_data_atmo(i,1) > 0
    i = i + 1;
end
autocorr_time_lag_atmo(i:end) = [];
autocorr_data_atmo(i:end) = [];

integral_length_atmo = trapz(autocorr_time_lag_atmo,autocorr_data_atmo)*r_atmo; 

%   dataCenter
timeLags_dataCenter = length(dataCenterline)-1;
autocorr_time_lag_dataCenter = (1-1:timeLags_dataCenter)';
autocorr_data_dataCenter = autocorr(dataCenterline,timeLags_dataCenter);
%delete data lower than 0 of the autocorrelation (visually dicided)
i = 1;
while autocorr_data_dataCenter(i,1) > 0
    i = i + 1;
end
autocorr_time_lag_dataCenter(i:end) = [];
autocorr_data_dataCenter(i:end) = [];

integral_length_dataCenter = trapz(autocorr_time_lag_dataCenter,autocorr_data_dataCenter)*r_dataCenter; 
%% kolmogorov length
% atmo
epsilon_atmo = nanmean(atmosphere).^3 / integral_length_atmo;
visco_air = 1.5*10^-5;
kolmogorov_length_atmo = (visco_air^3/epsilon_atmo)^(1/4);
% dataCenter
epsilon_dataCenter = nanmean(dataCenterline).^3 / integral_length_dataCenter;
visco_air = 1.5*10^-5;
kolmogorov_length_dataCenter = (visco_air^3/epsilon_dataCenter)^(1/4);
%% Velocity increments for r = 2^m for m = 1-8
lag_array = [2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8];
% atmo
lag_increments_atmo = lagmatrix(atmosphere(:,1),[0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8]);

count = 1;
for i = 2:length(lag_increments_atmo(1,:))
    disp(i);
    increments_atmo(:,count) = lag_increments_atmo(:,i) - lag_increments_atmo(:,1);
    count = count + 1
end
increments_atmo_space = increments_atmo * r_atmo;
% dataCenter
lag_increments_center = lagmatrix(dataCenterline(:,1),[0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8]);

count = 1;
for i = 2:length(lag_increments_center(1,:))
    disp(i);
    increments_center(:,count) = lag_increments_center(:,i) - lag_increments_center(:,1);
    count = count + 1
end
increments_center_space = increments_center * r_dataCenter;
%% determine structure function <u_r^2>
% atmo
structure_function_square_atmo = nanmean(increments_atmo.^2);
% dataCenter
structure_function_square_center = nanmean(increments_center.^2);
%% determine strucutre function <u_r^n>
n = [1,2,3,4,5,6]
% atmo
for i = n 
    structure_function_atmo(:,i) = nanmean(increments_atmo.^i);
end
% dataCenter
for i = n 
    structure_function_center(:,i) = nanmean(increments_center.^i);
end  
%% Plot of structure functions
% atmo
ylim_atmo = max(structure_function_atmo(:))
x_pos = integral_length_atmo/r_atmo;
figure
%plot([x_pos x_pos], [0 ylim]);
hold on
for i = 1:length(structure_function_atmo(1,:))
    loglog(lag_array,abs(structure_function_atmo(:,i)),'-s');
    grid on
end
hold off
title('Structure functions  of atmosphere')
xlabel('Lags [m]')
ylabel('Stucturefunctions u(r,n)')
% dataCenter
figure
ylim_center = max((abs(structure_function_center(:))))
x_pos_center = integral_length_dataCenter;
for i = 1:length(structure_function_center(1,:))
    loglog(lag_array,abs(structure_function_center(:,i)),'-s');
    hold on
    grid on
end
loglog([x_pos_center x_pos_center], [0.000001 ylim_center]);
ylim([0.000001,ylim_center])
hold off
title('Structure functions  of dataCenter')
xlabel('Lags [m]')
ylabel('Stucturefunctions u(r,n)')

%% PDF of increments
%   atmo
for i = 1:length(increments_atmo_space(1,:))
    [hist_y, hist_x] = hist(increments_atmo_space(:,i),50);
    hist_y = hist_y/sum(hist_y);
    hist_y_array(:,i) = hist_y;
    hist_x_array(:,i) = hist_x;   
end
figure
for i = 1:length(increments_atmo_space(1,:))
    hold on
    semilogy(hist_x_array(:,i),hist_y_array(:,i))  
    set(gca,'yscale','log')
end
hold off
title('PDF u increment of atmosphere')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
% dataCenter
for i = 1:length(increments_center_space(1,:))
    [hist_y, hist_x] = hist(increments_center_space(:,i),50);
    hist_y = hist_y/sum(hist_y);
    hist_y_array_center(:,i) = hist_y;
    hist_x_array_center(:,i) = hist_x;   
end
figure
for i = 1:length(increments_center_space(1,:))
    hold on
    semilogy(hist_x_array_center(:,i),hist_y_array_center(:,i))  
    set(gca,'yscale','log')
end
hold off
title('PDF u increment of dataCenter')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
%% Estimate scaling exponents
