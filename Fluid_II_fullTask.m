%% 1 basic characteristic
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
disp('Computing 10min means and stddev');  
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
    end

    for i = 1:length(means_interval600_center(~isnan(means_interval600_center(:,j))))
        fluc_center((i-1)*timeInterval+1:i*timeInterval,k) = dataCenterline((i-1)*timeInterval+1:i*timeInterval,1) - means_interval600_center(i,j);
    end
    j = j + 2
    k = k +1
end
save('fluctuations.mat','fluc_atmo', 'fluc_center');
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

Tabele_mean= table(mean_60_atmo,mean_600_atmo,mean_6000_atmo,mean_comp_atmo,'RowNames',lines_name_tab);

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

Tabele_std= table(std_60_atmo,std_600_atmo,std_6000_atmo,std_comp_atmo,'RowNames',lines_name_tab);

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

Tabele_deg_turb = table(deg_turb_60_atmo,deg_turb_600_atmo,deg_turb_6000_atmo,deg_turb_comp_atmo,'RowNames',lines_name_tab);
%% Prob Densities
close all;
histo_fluc_std(1:timeInterval) = fluc_atmo(1:timeInterval)/means_interval600_atmo(1,2)

[fluc_plot,x] = hist(fluc_atmo(1:timeInterval),20);
fluc_plot = fluc_plot / sum(fluc_plot);
fit_fluc = fitdist(fluc_plot','Normal')
fit_pdf = pdf(fit_fluc,x);
bar(x,fluc_plot);
hold on;
plot(x,fit_pdf);
hold off;

%% Plotting

figure()
plot(atmosphere, '-r');
hold on;
for i = 1:length(means_interval600_atmo)
    plot_data_mean(1:timeInterval,1) = means_interval600_atmo(i,1);
    plot_data_mean(1:timeInterval,2) = means_interval600_atmo(i,1)+5*means_interval600_atmo(i,2);
    plot_data_mean(1:timeInterval,3) = means_interval600_atmo(i,1)-5*means_interval600_atmo(i,2);
    plot(((i-1)*timeInterval+1:i*timeInterval), plot_data_mean(1:timeInterval,1),'-g') 
    plot(((i-1)*timeInterval+1:i*timeInterval), plot_data_mean(1:timeInterval,2),'-b') 
    plot(((i-1)*timeInterval+1:i*timeInterval), plot_data_mean(1:timeInterval,3),'-b') 
end
hold off

