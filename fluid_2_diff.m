%Import data
%Lab data
data = dlmread('jfm_data_block1.txt','\t');
data_name = 'Labordaten';
fq = 8000; %sample frequencies
axis_5_3_PSD = 10.3;
stopp_autocorr = 900;

%Real Data
%data = dlmread('AtmosphericData_July_fs10Hz_Kurz.txt','\t');
%data_name = 'Daten von Fino1';
%fq = 10; %sample frequencies
%axis_5_3_PSD = 11.0;
%stopp_autocorr = 950;
%nan_array = isnan(data);
%mean = nanmean(data);
%data(nan_array) = mean;

T = 1/fq;
ny_air = 1.42*10^-5;

lines = size(data,1); 

%power spectral density
Fs = fq;            % Sampling frequency  
L = lines*T*1000;   % Length of signal in ms 

Y = fft(data);

P2 = abs(Y);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
PSD = 2*P1.^2;

f = Fs*(0:(L/2))/L;

%figure
%plot(f,PSD)
%title('One sided Amplitude Spectrum of u(t)')
%xlabel('f (Hz)')
%ylabel('|P1(f)|')

P_gerade = 10.^(-5/3 *log10(f)+axis_5_3_PSD);
yy = smooth(PSD,20); 

figure
loglog(f,PSD)        %regular spectra
hold on
loglog(f,yy)         %smootrhed spectra
plot(f,P_gerade,'Color','black')     %scaling   
title('Spektrale Leistungsdichte')
xlabel('f (Hz)')
ylabel('P(f)')
xlim([0 size(f,2)])
legend('Reguläre PSD','Geglättete PSD','k^-^5^/^3 Skalierung');
saveas (gcf,['Spektrale Leistungsdichte, ',data_name,'.png']);


%integral length/autocorr
timeLags = 125000-1;
autocorr_time_lag = (1-1:timeLags)';
autocorr_data = autocorr(data,timeLags);
%delete data lower than 0 of the autocorrelation (visually dicided)
autocorr_data_clean = autocorr_data;
autocorr_time_lag_clean=autocorr_time_lag;
autocorr_time_lag_clean(stopp_autocorr:end) = [];
autocorr_data_clean(stopp_autocorr:end) = [];
figure
hold on
plot(autocorr_time_lag,autocorr_data);
xlabel('Zeitschritt')
ylabel('Korrelation')
title('Autokorrelation, gesamter Bereich')
hold off

figure
plot(autocorr_time_lag_clean,autocorr_data_clean);
xlabel('Zeitschritt')
ylabel('Korrelation')
title('Autokorrelation, relevanter Bereich')
saveas (gcf,['Autokorrelation, ',data_name,'.png']);

%JOINT PROBABILITY


%Correlated
time_lag_choosen = 6;
u_of_t = data(1:end-time_lag_choosen)-nanmean(data);
u_of_t_plus_tau =data(time_lag_choosen+1:end)-nanmean(data);
deci = 1; %decimals
u_round = round(u_of_t,deci);
u_t_round = round(u_of_t_plus_tau,deci);

% max(u_round) min(u_round)
u_xarray = min(u_round):10^-deci: max(u_round);
u_yarray = min(u_t_round):10^-deci: max(u_t_round);
edgesy = u_yarray - 10^-deci/2;
edgesy(length(edgesy)+1)=edgesy(length(edgesy))+10^-deci/2;

figure
h = histogram(u_of_t,edgesy);
counts_u_t = h.Values;
prob_vector_u_t = counts_u_t./sum(counts_u_t);

h = histogram(u_of_t_plus_tau,edgesy);
counts_u_t_tau = h.Values;
prob_vector_u_t_tau = counts_u_t_tau./sum(counts_u_t_tau);

close Figure 4

for i = 1:length(prob_vector_u_t)
    
help_array = prob_vector_u_t(i)*prob_vector_u_t_tau;    
prob_field(:,i)=help_array';
end


figure
contourf(u_xarray,u_yarray,prob_field);
title(['Verbundwahrscheinlichkeitsdichte, Daten korreliert, ',data_name], 'FontSize', 15)
xlabel('u\prime(t) [m/s]', 'FontSize', 15)
ylabel('u\prime(t+\tau)[m/s]', 'FontSize', 15)
colorbar('Ticks',[0,0.001,0.005,0.01], 'TickLabels',{'0','0.001','0.005','0.01'})
saveas (gcf,['Joint_prob_correlated u, ',data_name,'.png']);








%NotCorrelated
time_lag_choosen = 2000;
u_of_t_2 = data(1:end-time_lag_choosen)-nanmean(data);
u_of_t_plus_tau_2 =data(time_lag_choosen+1:end)-nanmean(data);
deci = 1; %decimals
u_round_2 = round(u_of_t_2,deci);
u_t_round_2 = round(u_of_t_plus_tau_2,deci);

% max(u_round) min(u_round)
u_xarray_2 = min(u_round_2):10^-deci: max(u_round_2);
u_yarray_2 = min(u_t_round_2):10^-deci: max(u_t_round_2);
edgesy_2 = u_yarray_2 - 10^-deci/2;
edgesy_2(length(edgesy_2)+1)=edgesy_2(length(edgesy_2))+10^-deci/2;

figure
h = histogram(u_of_t_2,edgesy_2);
counts_u_t_2 = h.Values;
prob_vector_u_t_2 = counts_u_t_2./sum(counts_u_t_2);

h = histogram(u_of_t_plus_tau_2,edgesy);
counts_u_t_tau_2 = h.Values;
prob_vector_u_t_tau_2 = counts_u_t_tau_2./sum(counts_u_t_tau_2);

for i = 1:length(prob_vector_u_t_2)
    
help_array = prob_vector_u_t_2(i)*prob_vector_u_t_tau_2;    
prob_field_2(:,i)=help_array';
end

close Figure 5

figure
contourf(u_yarray_2,u_yarray_2,prob_field_2);
title(['Verbundwahrscheinlichkeitsdichte, Daten nicht korreliert, ',data_name], 'FontSize', 15)
xlabel('u\prime(t) [m/s]', 'FontSize', 15)
ylabel('u\prime(t+\tau)[m/s]', 'FontSize', 15)
colorbar('Ticks',[0,0.001,0.005,0.01], 'TickLabels',{'0','0.001','0.005','0.01'})
saveas (gcf,['Joint_prob_not_correlated u, ',data_name,'.png']);





%INTEGRAL LENGTH/KOLMOGOROV LENGTH

taylor_hypothesis_factor = 1/fq * nanmean(data); %frozen turbulence

%Integral length
integral_length = trapz(autocorr_time_lag(1:stopp_autocorr),autocorr_data(1:stopp_autocorr))*taylor_hypothesis_factor;

%Kolmogorov
u_of_x = data(1:end-1);
u_of_xplustau = data(1+1:end);
increment_kolmo = u_of_x - u_of_xplustau ;

sec_derevative = increment_kolmo-nanmean(data);
sec_derevative = (sec_derevative./taylor_hypothesis_factor).^2; 
sec_der_av = nanmean(sec_derevative);
eta_kolmo = (ny_air^2/(15*sec_der_av))^(1/4);


%CALCULATE INCREMENTS
m_max = 10;
m = 0: m_max;
increment_lags = 2.^m;

increment_array = nan (length(data),m_max+1);
for i = 1 : m_max+1
    
u_of_x = data(1:end-increment_lags(i));
u_of_xplustau = data(1+increment_lags(i):end);

increment_array(1:length(u_of_x),i) = u_of_x - u_of_xplustau ;

end
clear u_of_x u_of_xplustau i

figure
hold on
histogram(increment_array(:,8),50,'Normalization','probability','FaceColor','g')
histogram(increment_array(:,4),50,'Normalization','probability','FaceColor','b')
histogram(increment_array(:,1),50,'Normalization','probability','FaceColor','k')
title(['Verteilung der Inkremente, ',data_name], 'FontSize', 15)
xlabel('ur [m/s]', 'FontSize', 15)
ylabel('Relative Häufigkeit [1]', 'FontSize', 15)
legend('m = 7','r  = 3','r  = 0')
hold off
saveas (gcf,['Inkremente, ',data_name,'.png']);


%STRUCTURE FUNCTION
orders_cons = 6;
for i = 2 :orders_cons    
incre_arr_mom = increment_array.^i;     
incre_arr_mom_mean =nanmean(incre_arr_mom);  
structure(i,:)= abs(incre_arr_mom_mean);
end


figure
hold on
for i = 2 :orders_cons
scatter(increment_lags,structure(i,:)) 
set(gca,'xscale','log','yscale','log')
end
hold off


structure_function_square = nanmean(increment_array.^2);
r_increments = increment_lags .* taylor_hypothesis_factor;
u_prime = data-nanmean(data);


taylor_length_function = nanmean(u_prime.^2)./structure_function_square .* r_increments .^2;

figure
scatter(r_increments,taylor_length_function)
title('Taylor length function')
xlabel('r [m]')
ylabel('lamda^2 [m^2]')


