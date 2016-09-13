%% two-point quantities
% load data
clear all;
close all;
atmosphere = load ('AtmosphericData_July_fs10Hz_Kurz.txt');
dataCenterline = load('jfm_data_block1.txt');
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

Fs_center = 8000% sampling frequency dataCenter

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
%% PSD Smoothing
PSD_center_smooth = smooth(PSD_center,20);
PSD_atmo_smooth = smooth(PSD_atmo,20);


%% fit for 5/3
P_gerade_atmo = 10.^(-5/3 *log10(f_atmo)+4);
P_gerade_atmo_fluc = 10.^(-5/3 *log10(f_atmo_fluc)+4);
P_gerade_center = 10.^(-5/3 *log10(f_center)+8.5);
P_gerade_center_fluc = 10.^(-5/3 *log10(f_center_fluc)+10.3);

% fit plot atmosphere u
figure
loglog(f_atmo',PSD_atmo);
hold on
plot(f_atmo',P_gerade_atmo)
loglog(f_atmo',PSD_atmo_smooth)
title('PSD atmo')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_atmo,1)])
saveas(gcf, 'report/figures/power_spectrum_atmo_smooth.png');

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
loglog(f_center',PSD_center_smooth);
plot(f_center',P_gerade_center)
title('PSD dataCenter')
xlabel('f (Hz)')
ylabel('Power spectral density')
xlim([0 size(f_center',1)])
saveas(gcf, 'report/figures/power_spectrum_center_smooth.png');


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
saveas(gcf,'report/figures/atmo_numerically_equal.png')
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
saveas(gcf,'report/figures/center_numerically_equal.png')

%%  Joint probability dist
% atmosphere , create lagmatrix
xlag_atmo = lagmatrix(fluc_atmosphere(:,1),[0 1 5 10 100]);
% check for correlation
xlag_names = [1, 5, 10, 100];

figure
for i = 2 : length(xlag_atmo(1,:))
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
    subplot(4,1,i-1)
    contour(xi,yi,pdfData);
    title(['shift of ' num2str(xlag_names(i-1))])
end
saveas(gcf,'report/figures/jpdf_atmo.png')
% dataCenterline , create lagmatrix
xlag_dataCenter = lagmatrix(fluc_dataCenter(:,1),[0 1 5 10 100]);
% check for correlation

figure
for i = 2 : length(xlag_dataCenter(1,:))
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
    subplot(4,1,i-1)
    contour(xi,yi,pdfData);
    title(['shift of ' num2str(xlag_names(i-1))])
end
saveas(gcf,'report/figures/jpdf_center.png')



%% Taylor's hypothesis

r_atmo = nanmean(atmosphere) * 1/Fs_atmo; % diameter of frozen turbulence structure
r_dataCenter = nanmean(dataCenterline) * 1/Fs_center; % diameter of frozen turbulence structure
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

figure
plot(autocorr_data_atmo)
xlabel('Lags [m]')
ylabel('Correlation Coefficient')
saveas(gcf,'corr_atmo_int.png')

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
figure
plot(autocorr_data_dataCenter)
xlabel('Lags [m]')
ylabel('Correlation Coefficient')
saveas(gcf,'corr_atmo_int.png')

integral_length_dataCenter = trapz(autocorr_time_lag_dataCenter,autocorr_data_dataCenter)*r_dataCenter; 
%% kolmogorov length
% atmo
ny_air = 1.42*10^-5;
% alternative
u_of_x = atmosphere(1:end-1);
u_of_xplustau = atmosphere(1+1:end);
increment_kolmo = u_of_x - u_of_xplustau ;

sec_derevative = increment_kolmo-nanmean(atmosphere);
sec_derevative = (sec_derevative./r_atmo).^2; 
sec_der_av_atmo = nanmean(sec_derevative);
eta_kolmo_atmo = (ny_air^2/(15*sec_der_av_atmo))^(1/4);

% lab
u_of_x = dataCenterline(1:end-1);
u_of_xplustau = dataCenterline(1+1:end);
increment_kolmo = u_of_x - u_of_xplustau ;

sec_derevative = increment_kolmo-nanmean(dataCenterline);
sec_derevative = (sec_derevative./r_dataCenter).^2; 
sec_der_av_center = nanmean(sec_derevative);
eta_kolmo_center = (ny_air^2/(15*sec_der_av_center))^(1/4);

%

epsilon_atmo = 15*ny_air*sec_der_av_atmo;
epsilon_center = 15*ny_air*sec_der_av_center;

%kolmogorov_length_atmo = (visco_air^3/epsilon_atmo)^(1/4);
% dataCenter
epsilon_dataCenter = nanmean(dataCenterline).^3 / integral_length_dataCenter;
visco_air = 1.42*10^-5;
%kolmogorov_length_dataCenter = (visco_air^3/epsilon_dataCenter)^(1/4);
%% Taylor Length
taylor_length_atmo = sqrt(nanmean(fluc_atmosphere.^2)/sec_der_av_atmo);
taylor_length_center = sqrt(nanmean(fluc_dataCenter.^2)/sec_der_av_center);
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
r_lag_atmo = lag_array * r_atmo;

FigH = figure
subplot(3,1,1)
plot(increments_atmo(1:2000,1))
xlabel('Time in [0.1 s]')
%ylabel('U_x fluctuations in [m/s]')
title('Spatial Series of Increments from atmosphere for r^1')

subplot(3,1,2)
plot(increments_atmo(1:2000,2))
xlabel('Time in [0.1 s]')
%ylabel('U_x fluctuations in [m/s]')
title('Spatial Series of Increments from atmosphere for r^2')

subplot(3,1,3)
plot(increments_atmo(1:2000,5))
xlabel('Time in [0.1 s]')
title('Spatial Series of Increments from atmosphere for r^5')

AxesH    = findobj(FigH, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'U_x in [m/s]')



saveas(gcf,'report/figures/increments_array_atmo.png')


% dataCenter
lag_increments_center = lagmatrix(dataCenterline(:,1),[0,2^1,2^2,2^3,2^4,2^5,2^6,2^7,2^8]);

count = 1;
for i = 2:length(lag_increments_center(1,:))
    disp(i);
    increments_center(:,count) = lag_increments_center(:,i) - lag_increments_center(:,1);
    count = count + 1
end

FigH = figure
subplot(3,1,1)
plot(increments_center(1:2000,1))

xlabel('Time in [0.125 ms]')
title('Spatial Series of Increments from JFM for r^1')

subplot(3,1,2)
plot(increments_center(1:2000,2))
xlabel('Time in [0.125 ms]')
title('Spatial Series of Increments from JFM for r^2')

subplot(3,1,3)
plot(increments_center(1:2000,5))
xlabel('Time in [0.125 ms]')
title('Spatial Series of Increments from JFM for r^5')

AxesH    = findobj(FigH, 'Type', 'Axes');
YLabelHC = get(AxesH, 'YLabel');
YLabelH  = [YLabelHC{:}];
set(YLabelH, 'String', 'U_x in [m/s]')
saveas(gcf,'report/figures/increments_array_center.png')

r_lag_center = lag_array * r_dataCenter;
%% determine structure function <u_r^2>
% atmo
structure_function_square_atmo = nanmean(increments_atmo.^2);
% dataCenter
structure_function_square_center = nanmean(increments_center.^2);
%% determine strucutre function <u_r^n>
n = [1,2,3,4,5,6,7,8]
% atmo
for i = n 
    structure_function_atmo(:,i) = abs(nanmean(increments_atmo.^i));
    structure_function_atmo_ess(:,i) = nanmean(abs(increments_atmo.^i));
end
% dataCenter
for i = n 
    structure_function_center(:,i) = abs(nanmean(increments_center.^i));
    structure_function_center_ess(:,i) = nanmean(abs(increments_center.^i));
end  
%% Taylor Length Function
taylor_length_function_atmo = nanmean(atmosphere.^2)./structure_function_square_atmo .*  r_lag_atmo.^2;
taylor_length_function_center = nanmean(dataCenterline.^2)./structure_function_square_center .*  r_lag_center.^2;

figure
scatter(r_lag_atmo,taylor_length_function_atmo)
title('Taylor length function')
xlabel('r [m]')
ylabel('lamda^2 [m^2]')
%% Plot of structure functions
% atmo
ylim_atmo = max(structure_function_atmo(:))
x_pos_tay_atmo = taylor_length_atmo*100;
figure
for i = 1:length(structure_function_atmo(1,:))
    loglog(r_lag_atmo,abs(structure_function_atmo(:,i)),'-s');
    hold on
    grid on
end
loglog([x_pos_tay_atmo x_pos_tay_atmo], [0.00001 ylim_atmo]);
hold off
text(x_pos_tay_atmo,ylim_atmo/2^12,'Taylor length')
title('Structure functions  of atmosphere')
xlabel('Lags [m]')
ylabel('Stucturefunctions u(r,n)')
legend('Structure Order 1','Structure Order 2','Structure Order 3','Structure Order 4','Structure Order 5','Structure Order 6','Location','northwest')
saveas(gcf,'report/figures/structure_atmo.png')


% dataCenter
figure
ylim_center = max((abs(structure_function_center(:))))
x_pos_center = integral_length_dataCenter;
x_pos_kolmo_center = eta_kolmo_center;
x_pos_tay_center = 0.00485;
for i = 1:length(structure_function_center(1,:))
    loglog(r_lag_center,abs(structure_function_center(:,i)),'-s');
    hold on
    grid on
end
loglog([x_pos_center x_pos_center], [0.000001 ylim_center]);
loglog([x_pos_tay_center x_pos_tay_center], [0.000001 ylim_center]);
hold off
title('Structure functions  of ')
xlabel('Lags [m]')
ylabel('Stucturefunctions u(r,n)')
text(x_pos_tay_center,ylim_center/2^12,'Taylor length')
text(x_pos_center,ylim_center/2^12,'Integral length')
legend('Structure Order 1','Structure Order 2','Structure Order 3','Structure Order 4','Structure Order 5','Structure Order 6','Location','northwest')
saveas(gcf,'report/figures/structure_center.png')

%% PDF of increments
%   atmo
for i = 1:length(increments_atmo(1,:))
    [hist_y, hist_x] = hist(increments_atmo(:,i),125);
    hist_y = hist_y/sum(hist_y);
    hist_y_array(:,i) = hist_y;
    hist_x_array(:,i) = hist_x;   
end
figure
for i = 1:length(increments_atmo(1,:))
    hold on
    semilogy(hist_x_array(:,i),hist_y_array(:,i))  
    set(gca,'yscale','log')
end
hold off
title('PDF u increment of atmosphere')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 0.1])
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7','m=8')
saveas(gcf,'report/figures/pdf_increments_atmo.png')
% dataCenter
for i = 1:length(increments_center(1,:))
    [hist_y, hist_x] = hist(increments_center(:,i),125);
    hist_y = hist_y/sum(hist_y);
    hist_y_array_center(:,i) = hist_y;
    hist_x_array_center(:,i) = hist_x;   
end
figure
for i = 1:length(increments_center(1,:))
    hold on
    semilogy(hist_x_array_center(:,i),hist_y_array_center(:,i))  
    set(gca,'yscale','log')
end
hold off
title('PDF u increment of dataCenter')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 0.1])
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7','m=8')
saveas(gcf,'report/figures/pdf_increments_center.png')

%   atmo
for i = 1:length(increments_atmo(1,:))
    [hist_y, hist_x] = hist(increments_atmo(:,i),125);
    hist_y = hist_y/sum(hist_y);
    hist_y = hist_y/nanstd(increments_atmo(:,i));
    hist_y_array(:,i) = hist_y;
    hist_x_array(:,i) = hist_x;   
end
figure
for i = 1:length(increments_atmo(1,:))
    hold on
    semilogy(hist_x_array(:,i),hist_y_array(:,i))  
    set(gca,'yscale','log')
end
hold off
title('Normed PDF u increment of atmosphere')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 1])
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7','m=8')
saveas(gcf,'report/figures/pdf_increments_atmo_normed.png')
% dataCenter
for i = 1:length(increments_center(1,:))
    [hist_y, hist_x] = hist(increments_center(:,i),125);
    hist_y = hist_y/sum(hist_y);
    hist_y = hist_y/nanstd(increments_center(:,i));
    hist_y_array_center(:,i) = hist_y;
    hist_x_array_center(:,i) = hist_x;   
end
figure
for i = 1:length(increments_center(1,:))
    hold on
    semilogy(hist_x_array_center(:,i),hist_y_array_center(:,i))  
    set(gca,'yscale','log')
end
hold off
title('Normed PDF u increment of dataCenter')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 10])
legend('m=1','m=2','m=3','m=4','m=5','m=6','m=7','m=8');
saveas(gcf,'report/figures/pdf_increments_center_normed.png')

% further investigation
figure
[hist_y, hist_x] = hist(increments_center(:,1),125);
hist_y = hist_y/sum(hist_y);
hist_y = hist_y/nanstd(increments_center(:,1));
hist_y_array_center(:,1) = hist_y;
hist_x_array_center(:,1) = hist_x;  
f = fit(hist_x_array_center(:,1),hist_y_array_center(:,1),'gauss2')
plot(f,hist_x_array_center(:,1),hist_y_array_center(:,1))
set(gca,'yscale','log')
title('Normed PDF u increment of dataCenter ')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 10])
legend('m=1');
saveas(gcf,'report/figures/pdf_increments_center_normed_1.png')

% further investigation
figure
[hist_y, hist_x] = hist(increments_center(:,8),125);
hist_y = hist_y/sum(hist_y);
hist_y = hist_y/nanstd(increments_center(:,8));
hist_y_array_center(:,8) = hist_y;
hist_x_array_center(:,8) = hist_x;  
f = fit(hist_x_array_center(:,8),hist_y_array_center(:,8),'gauss2')
plot(f,hist_x_array_center(:,8),hist_y_array_center(:,8))
set(gca,'yscale','log')
title('Normed PDF u increment of dataCenter ')
xlabel('Fluctuation U_r [m/s]')
ylabel('Probability Density P(U_r) [1]')
ylim([0.001 10])
legend('m=1');
saveas(gcf,'report/figures/pdf_increments_center_normed_8.png')
%% Estimate scaling exponents
% at LAG 3 % LN BENUTZEN!!!
for i = 1:length(structure_function_center(1,:))
    slope_center(i) = (log10(structure_function_center(end,i))-log10(structure_function_center(4,i)))/(log10(r_lag_center(end))-log10(r_lag_center(4)));
    slope_center_ess(i) = (log10(structure_function_center_ess(end,i))-log10(structure_function_center_ess(4,i)))/(log10(r_lag_center(end))-log10(r_lag_center(4)));
    kolmo41_slope_center(i) = i/3*log10(epsilon_center);
    kolmo62_slope_center(i) = i/3 - 0.004*(i*(i-3))/18
end

for i = 1:length(structure_function_atmo(1,:))
    slope_atmo(i) = (log10(structure_function_atmo(end,i))-log10(structure_function_atmo(4,i)))/(log10(r_lag_atmo(end))-log10(r_lag_atmo(4)));
    slope_atmo_ess(i) = (log10(structure_function_atmo_ess(end,i))-log10(structure_function_atmo_ess(4,i)))/(log10(r_lag_atmo(end))-log10(r_lag_atmo(4)));
    kolmo41_slope_atmo(i) = i/3*abs(log10(epsilon_atmo));
    kolmo62_slope_atmo(i) = i/3 - 0.004*(i*(i-3))/18
end

figure
plot(slope_center)
hold on
plot(slope_center_ess)
plot(kolmo41_slope_center)
plot(kolmo62_slope_center)
hold off
title('Scaling Exponents Function of dataSet')
xlabel('Exponent n')
ylabel('Scaling Exponents')
legend('Scaling Exponents','Scaling Exponents (ESS)', 'Kol41 Scaling', 'Kol62 Scaling')
saveas(gcf,'report/figures/scale_compare_center.png')

figure
plot(slope_atmo)
hold on
plot(slope_atmo_ess)
plot(kolmo41_slope_atmo)
plot(kolmo62_slope_atmo)
hold off
title('Scaling Exponents Function of atmosphere')
xlabel('Exponent n')
ylabel('Scaling Exponents')
legend('Scaling Exponents','Scaling Exponents (ESS)', 'Kol41 Scaling', 'Kol62 Scaling')
saveas(gcf,'report/figures/scale_compare_atmo.png')
%% ESS
figure
for i = 1:length(structure_function_center(1,:))
    loglog(structure_function_center_ess(:,3),structure_function_center_ess(:,i))
    hold on
end
xlabel('Structure Function S(3,r)')
ylabel('Strucuture Function S(n,r)')
saveas(gcf,'report/figures/ess_compare_center.png')
hold off

figure
for i = 1:length(structure_function_atmo(1,:))
    loglog(structure_function_atmo_ess(:,3),structure_function_atmo_ess(:,i))
    hold on
end
xlabel('Structure Function S(3,r)')
ylabel('Strucuture Function S(n,r)')
saveas(gcf,'report/figures/ess_compare_atmo.png')
hold off

%% n - point quantities
xlag_atmo = increments_atmo;
xlag_center = increments_center;

figure
for i = 1 : 6
    vector_n(:,1) = xlag_atmo(:,1);
    vector_n(:,2) = xlag_atmo(:,i);
    xi = linspace(min(vector_n(:,1)), max(vector_n(:,1)), 200);
    yi = linspace(min(vector_n(:,2)), max(vector_n(:,2)), 200);
    hst = hist3(vector_n(:,1:2),{xi yi}); %removed extra '
    hst_sim = hist(vector_n(:,2), 200);
    hst_sim = meshgrid(hst_sim,hst_sim);
    %normalize
    dx = xi(2)-xi(1);
    dy = yi(2)-yi(1);
    area = dx*dy;
    pdfData = hst/sum(sum(hst))/area;
    pdfData_new = pdfData/hst_sim;
    subplot(3,2,i)
    s1 = contour(xi,yi,pdfData);
    xlabel('Increment u_1');
    ylabel('Increment u_2');
    legend(['r=',num2str(r_atmo*2.^i)],'location','northwest')
end
saveas(gcf,'report/figures/npoint_atmo.png')

figure
for i = 1 : 6
    vector_new(:,1) = xlag_center(:,1);
    vector_new(:,2) = xlag_center(:,i);
    xi = linspace(min(vector_new(:,1)), max(vector_new(:,1)), 200);
    yi = linspace(min(vector_new(:,2)), max(vector_new(:,2)), 200);
    hst = hist3(vector_new(:,1:2),{xi yi}); %removed extra '
    hst_sim = hist(vector_new(:,2), 200);
    %hst_sim = meshgrid(hst_sim,hst_sim);
    %normalize
    dx = xi(2)-xi(1);
    dy = yi(2)-yi(1);
    area = dx*dy;
    pdfData = hst/sum(sum(hst))/area;
    pdfData_new = pdfData/hst_sim;
    subplot(3,2,i)
    contour(xi,yi,pdfData);
    xlabel('Increment u_1');
    ylabel('Increment u_2');
    legend(['r=',num2str(r_dataCenter*2.^i)],'location','northwest')
end
saveas(gcf,'report/figures/npoint_center.png')