clear all;
close all;

%Import data
data_lab = dlmread('jfm_data_block1.txt','\t');

%sample frequencies
fq_lab = 8000;
ny_air = 1.42*10^-5;

lines = size(data_lab,1);

%power spectral density

Fs = fq_lab;            % Sampling frequency  
L = lines;             % Length of signal

Y = fft(data_lab);

P2 = abs(Y);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
PSD = 2*P1.^2;

f = Fs*(0:(L/2))/L;
figure
plot(f,PSD)
title('One sided Amplitude Spectrum of u(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
P_gerade = 10.^(-5/3 *log10(f)+8.8);



figure
loglog(f,PSD)
hold on
plot(f,P_gerade)
title('PSD')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f,2)])


yy = smooth(PSD,20);
figure
loglog(f,yy)
hold on
plot(f,P_gerade)
title('smoothed PSD')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f,2)])


%integral length
timeLags = 12500;
autocorr_time_lag = (1-1:timeLags)';
autocorr_data = autocorr(data_lab,timeLags);
%delete data lower than 0 of the autocorrelation (visually dicided)
autocorr_time_lag(950:end) = [];
autocorr_data(950:end) = [];

taylor_hypothesis_factor = 1/fq_lab * nanmean(data_lab); %frozen turbulence; diameter of structure

integral_length = trapz(autocorr_time_lag,autocorr_data)*taylor_hypothesis_factor;
eta_kolmo = ((ny_air^3 * integral_length)/(nanmean(data_lab)^3))^(1/4); 

%calculate increments
m_max = 8;
m = 0: m_max;
increment_lags = 2.^m;

increment_array = nan (length(data_lab),9);
for i = 1 : m_max+1
    
u_of_x = data_lab(1:end-increment_lags(i));
u_of_xplustau = data_lab(1+increment_lags(i):end);

increment_array(1:length(u_of_x),i) = u_of_x - u_of_xplustau ;

end
clear u_of_x u_of_xplustau i

structure_function_square = nanmean(increment_array.^2);

r_increments = increment_lags .* taylor_hypothesis_factor;
u_prime = data_lab-nanmean(data_lab);


taylor_length_function = nanmean(u_prime.^2)./structure_function_square .* r_increments .^2;

figure
scatter(r_increments,taylor_length_function)
title('Taylor length function')
xlabel('r [m]')
ylabel('lamda^2 [m^2]')


% m = moment(X,order)