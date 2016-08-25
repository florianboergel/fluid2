%% 1 basic characteristic
close all;
clear all;
%% 1.1 load data
atmosphere = load ('AtmosphericData_July_fs10Hz_Kurz.txt');
dataCenterline = load('Data_Centerline_FractalGrid_fs60kHz.txt');
%%
timeInterval = 60;

% mean value
mean_atmosphere = nanmean(atmosphere);
mean_dataCenterline = nanmean(dataCenterline);
% fluctuation 
std_atmosphere = nanmean(atmosphere);
std_dataCenterline = nanstd(dataCenterline);
% degree of turbulence
degree_turb_atmosphere = std_atmosphere.^2/mean_atmosphere.^2;
degree_turb_dataCenter = std_dataCenterline.^2/mean_dataCenterline.^2;


% init matrix
disp('Computing different means ...');  
means_interval600_atmo = NaN((floor(length(atmosphere(:,1))/timeInterval)),6);
means_interval600_center = NaN((floor(length(dataCenterline(:,1))/timeInterval)),6);

% different mean intervalls
j = 1
k = 1
for timeInterval = [60 600 6000]
    disp(timeInterval);
    for i = 1:floor(length(atmosphere(:,1))/timeInterval)
       velocities = atmosphere((i-1)*timeInterval+1:i*timeInterval,1);
       mean = nanmean(velocities);
       std = nanstd(velocities);
       means_interval600_atmo(i,j) = mean;
       means_interval600_atmo(i,j+1) = std;
       deg_turb_atmo(i,k) = means_interval600_atmo(i,j+1).^2 / (means_interval600_atmo(i,j+1)).^2;
    end

    for i = 1:floor(length(dataCenterline(:,1))/timeInterval)
       %treat wind vanes
       velocities = dataCenterline((i-1)*timeInterval+1:i*timeInterval,1);
       mean = nanmean(velocities);
       std = nanstd(velocities);
       means_interval600_center(i,j) = mean;
       means_interval600_center(i,j+1) = std;
       deg_turb_center(i,k) = means_interval600_center(i,j+1).^2 / (means_interval600_center(i,j+1)).^2;
    end
    % calculate fluctuations 
    length(means_interval600_atmo(~isnan(means_interval600_atmo(:,j))))
    for i = 1:length(means_interval600_atmo(~isnan(means_interval600_atmo(:,j))))
        fluc_atmo((i-1)*timeInterval+1:i*timeInterval,k) = atmosphere((i-1)*timeInterval+1:i*timeInterval,1) - means_interval600_atmo(i,j);
        fluc_atmo_norm((i-1)*timeInterval+1:i*timeInterval,k) = fluc_atmo((i-1)*timeInterval+1:i*timeInterval,k)./means_interval600_atmo(i,j+1);
    end

    for i = 1:length(means_interval600_center(~isnan(means_interval600_center(:,j))))
        fluc_center((i-1)*timeInterval+1:i*timeInterval,k) = dataCenterline((i-1)*timeInterval+1:i*timeInterval,1) - means_interval600_center(i,j);
        fluc_center_norm((i-1)*timeInterval+1:i*timeInterval,k) = fluc_center((i-1)*timeInterval+1:i*timeInterval,k)./means_interval600_center(i,j+1);
    end
    j = j + 2
    k = k +1
end

% u_dash
atmosphere_dash = atmosphere - mean_atmosphere;
dataCenterline_dash = dataCenterline - mean_dataCenterline;
save('fluctuations.mat','fluc_atmo', 'fluc_center','atmosphere_dash','dataCenterline_dash');
%% Comparison

%table: mean
lines_name_tab = {'min','max'};
mean_60_atmo = [min(means_interval600_atmo(:,1));max(means_interval600_atmo(:,1))];
mean_600_atmo = [min(means_interval600_atmo(:,3));max(means_interval600_atmo(:,3))];
mean_6000_atmo = [min(means_interval600_atmo(:,5));max(means_interval600_atmo(:,5))];
mean_comp_atmo = [mean_atmosphere;mean_atmosphere];
mean_60_center = [min(means_interval600_center(:,1));max(means_interval600_center(:,1))];
mean_600_center = [min(means_interval600_center(:,3));max(means_interval600_center(:,3))];
mean_6000_center = [min(means_interval600_center(:,5));max(means_interval600_center(:,5))];
mean_comp_center = [mean_dataCenterline;mean_dataCenterline];

Table_mean= table(mean_60_atmo,mean_600_atmo,mean_6000_atmo,mean_comp_atmo,mean_60_center,mean_600_center,mean_6000_center,mean_comp_center,'RowNames',lines_name_tab);

% table: std
lines_name_tab = {'min','max'};
std_60_atmo = [min(means_interval600_atmo(:,2));max(means_interval600_atmo(:,2))];
std_600_atmo = [min(means_interval600_atmo(:,4));max(means_interval600_atmo(:,4))];
std_6000_atmo = [min(means_interval600_atmo(:,6));max(means_interval600_atmo(:,6))];
std_comp_atmo = [std_atmosphere;std_atmosphere];
std_60_center = [min(means_interval600_center(:,2));max(means_interval600_center(:,2))];
std_600_center = [min(means_interval600_center(:,4));max(means_interval600_center(:,4))];
std_6000_center = [min(means_interval600_center(:,6));max(means_interval600_center(:,6))];
std_comp_center = [std_dataCenterline;std_dataCenterline];

Table_std= table(std_60_atmo,std_600_atmo,std_6000_atmo,std_comp_atmo,std_60_center,std_600_center,std_6000_center,std_comp_center,'RowNames',lines_name_tab);

%table: Degree of turbulence
lines_name_tab = {'min','max'};
deg_turb_60_atmo = [min(deg_turb_atmo(:,1));max(deg_turb_atmo(:,1))];
deg_turb_600_atmo = [min(deg_turb_atmo(:,2));max(deg_turb_atmo(:,2))];
deg_turb_6000_atmo = [min(deg_turb_atmo(:,3));max(deg_turb_atmo(:,3))];
deg_turb_comp_atmo = [degree_turb_atmosphere;degree_turb_atmosphere];
deg_turb_60_center = [min(deg_turb_center(:,1));max(deg_turb_center(:,1))];
deg_turb_600_center = [min(deg_turb_center(:,2));max(deg_turb_center(:,2))];
deg_turb_6000_center = [min(deg_turb_center(:,3));max(deg_turb_center(:,3))];
deg_turb_comp_center = [degree_turb_dataCenter;degree_turb_dataCenter];

Table_deg_turb = table(deg_turb_60_atmo,deg_turb_600_atmo,deg_turb_6000_atmo,deg_turb_comp_atmo,deg_turb_60_center,deg_turb_600_center,deg_turb_6000_center,deg_turb_comp_center,'RowNames',lines_name_tab);
%% statistics
close all;
% PDF u, full data
[f_u_tot_atmo,xi_u_tot_atmo] =  ksdensity (atmosphere);     
figure
plot (xi_u_tot_atmo,f_u_tot_atmo)
title('PDF, data from atmosphere.txt')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

[f_u_tot_center,xi_u_tot_center] =  ksdensity (dataCenterline);     
figure
plot (xi_u_tot_center,f_u_tot_center)
title('PDF, data from data_Centerline')
xlabel('wind speed [m/s]')
ylabel('probability [1]')


% PDF u', different means
[f_u_atmo_60,xi_u_60] =  ksdensity (fluc_atmo(:,1));     
figure
plot (xi_u_60,f_u_atmo_60)
title('PDF atmosphere, u_dash with timeInterval = 60')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

% PDF u', different means
[f_u_atmo_600,xi_u_600] =  ksdensity (fluc_atmo(:,2));     
figure
plot (xi_u_600,f_u_atmo_600)
title('PDF atmosphere, u_dash with timeInterval = 600')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

% PDF u', different means
[f_u_atmo_6000,xi_u_6000] =  ksdensity (fluc_atmo(:,3));     
figure
plot (xi_u_6000,f_u_atmo_6000)
title('PDF atmosphere, u_dash with timeInterval = 6000')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

% PDF u', different means
[f_u_center_60,xi_u_center_60] =  ksdensity (fluc_center(:,1));     
figure
plot (xi_u_center_60,f_u_center_60)
title('PDF dataCenter, u_dash with timeInterval = 60')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

% PDF u', different means
[f_u_center_600,xi_u_center_600] =  ksdensity (fluc_center(:,2));     
figure
plot (xi_u_center_600,f_u_center_600)
title('PDF dataCenter, u_dash with timeInterval = 600')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

% PDF u', different means
[f_u_center_6000,xi_u_center_6000] =  ksdensity (fluc_center(:,3));     
figure
plot (xi_u_center_6000,f_u_center_6000)
title('PDF dataCenter, u_dash with timeInterval = 6000')
xlabel('wind speed [m/s]')
ylabel('probability [1]')

%% Compare Plots
% u_dash
%   atmosphere
figure
plot (xi_u_60,f_u_atmo_60)
hold on
plot (xi_u_600,f_u_atmo_600)
plot (xi_u_6000,f_u_atmo_6000)
title('PDF atmosphere, u_dash comparison')
xlabel('wind speed [m/s]')
ylabel('probability [1]')
legend('averaging intervall = 60','averaging intervall = 600','averaging intervall = 6000')


%    dataCenter
figure
plot (xi_u_center_60,f_u_center_60)
hold on
plot (xi_u_center_600,f_u_center_600)
plot (xi_u_center_6000,f_u_center_6000)
title('PDF dataCenter, u_dash comparison')
xlabel('wind speed [m/s]')
ylabel('probability [1]')
legend('averaging intervall = 60','averaging intervall = 600','averaging intervall = 6000')

% u_dash normed by std
%   atmosphere
[f_u_norm_atmosphere_60,xi_u_norm_atmosphere_60] =  ksdensity (fluc_atmo_norm(:,1));     
[f_u_norm_atmosphere_600,xi_u_norm_atmosphere_600] =  ksdensity (fluc_atmo_norm(:,2));     
[f_u_norm_atmosphere_6000,xi_u_norm_atmosphere_6000] =  ksdensity (fluc_atmo_norm(:,3));  

figure
hold on
plot (xi_u_norm_atmosphere_60,f_u_norm_atmosphere_60)
plot (xi_u_norm_atmosphere_600,f_u_norm_atmosphere_600)
plot (xi_u_norm_atmosphere_6000,f_u_norm_atmosphere_6000)
hold off

title('PDF atmosphere - fluctuational part of velocity, data from laboratoty')
xlabel('wind speed fluctuation, normalized by standart deviation [1]')
ylabel('probability [1]')
legend('averaging intervall = 6000','averaging intervall = 600','averaging intervall = 60')


%   center
[f_u_norm_center_60,xi_u_norm_center_60] =  ksdensity (fluc_center_norm(:,1));     
[f_u_norm_center_600,xi_u_norm_center_600] =  ksdensity (fluc_center_norm(:,2));     
[f_u_norm_center_6000,xi_u_norm_center_6000] =  ksdensity (fluc_center_norm(:,3));     

figure
hold on
plot (xi_u_norm_center_60,f_u_norm_center_60)
plot (xi_u_norm_center_600,f_u_norm_center_600)
plot (xi_u_norm_center_6000,f_u_norm_center_6000)
hold off

title('PDF dataCenter - fluctuational part of velocity, data from laboratoty')
xlabel('wind speed fluctuation, normalized by standart deviation [1]')
ylabel('probability [1]')
legend('averaging intervall = 6000','averaging intervall = 600','averaging intervall = 60')


