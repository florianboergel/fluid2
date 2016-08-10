atmosphere = load ('AtmosphericData_July_fs10Hz_Kurz.txt');
dataCenterline = load('Data_Centerline_FractalGrid_fs60kHz.txt')

%% Basic characterisitic:
% mean value
mean_atmosphere = nanmean(atmosphere);
mean_dataCenterline = nanmean(dataCenterline);
% fluctuation 
std_atmosphere = nanmean(atmosphere);
std_dataCenterline = nanstd(dataCenterline);

% means and 10 min means
disp('Computing 10min means and stddev');  
means_interval10_atmo = NaN((floor(length(atmosphere(:,1))/600)),4);
means_interval10_center = NaN((floor(length(dataCenterline(:,1))/600)),4);

% Including degree of turbulence
for i = 1:floor(length(atmosphere(:,1))/600)
   %treat wind vanes
   velocities = atmosphere((i-1)*600+1:i*600,1);
   mean = nanmean(velocities);
   std = nanstd(velocities);
   means_interval10_atmo(i,1) = mean;
   means_interval10_atmo(i,2) = std;
   deg_turb_atmo(i,1) = means_interval10_atmo(i,2) / (means_interval10_atmo(i,1))
end

for i = 1:floor(length(dataCenterline(:,1))/600)
   %treat wind vanes
   velocities = dataCenterline((i-1)*600+1:i*600,1);
   mean = nanmean(velocities);
   std = nanstd(velocities);
   means_interval10_center(i,1) = mean;
   means_interval10_center(i,2) = std;
   deg_turb_center(i,1) = means_interval10_center(i,2) / (means_interval10_center(i,1))
end

% calculate fluctuations
for i = 1:length(means_interval10_atmo(:,1))
    fluc_atmo((i-1)*600+1:i*600,1) = atmosphere((i-1)*600+1:i*600,1) - means_interval10_atmo(i,1);
end

%% Prob Densities
close all;
histo_fluc_std(1:600) = fluc_atmo(1:600)/means_interval10_atmo(1,2)

[fluc_plot,x] = hist(fluc_atmo(1:600),20);
fluc_plot = fluc_plot / sum(fluc_plot);
fit_fluc = fitdist(fluc_plot','Normal')
fit_pdf = pdf(fit_fluc,x);
bar(x,fluc_plot);
hold on;
plot(x,fit_pdf);
hold off;

% %#METHOD 1: DIVIDE BY SUM
% figure()
% bar(x,f/trapz(x,f));hold on
% plot(x,g,'r');hold off
% title('FLUC')
% 
% 
% [f,x]=hist(histo_fluc_std(1:600),20);%# create histogram from a normal distribution.
% g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
% 
% %#METHOD 1: DIVIDE BY SUM
% figure()
% bar(x,f/trapz(x,f));hold on
% plot(x,g,'r');hold off
% title('FLUC/STD')
% 
% [f,x]=hist(atmosphere(1:600),20);%# create histogram from a normal distribution.
% g=1/sqrt(2*pi)*exp(-0.5*x.^2);%# pdf of the normal distribution
% 
% %#METHOD 1: DIVIDE BY SUM
% figure()
% bar(x,f/trapz(x,f));hold on
% plot(x,g,'r');hold off
% title('Velocity')

%% two-point quantities:
x = fluc_atmo;
Fs = 10; % Sampling frequency
t = 0:1/Fs:length(x)/10
fnyquist=Fs/2;N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

figure()
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
%% Plotting

figure()
plot(atmosphere, '-r');
hold on;
for i = 1:length(means_interval10_atmo)
    plot_data_mean(1:600,1) = means_interval10_atmo(i,1);
    plot_data_mean(1:600,2) = means_interval10_atmo(i,1)+5*means_interval10_atmo(i,2);
    plot_data_mean(1:600,3) = means_interval10_atmo(i,1)-5*means_interval10_atmo(i,2);
    plot(((i-1)*600+1:i*600), plot_data_mean(1:600,1),'-g') 
    plot(((i-1)*600+1:i*600), plot_data_mean(1:600,2),'-b') 
    plot(((i-1)*600+1:i*600), plot_data_mean(1:600,3),'-b') 
end
hold off
