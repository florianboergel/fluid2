clear all
close all

%Import data
%data = dlmread('jfm_data_block1.txt','\t');
%data_name = 'Labordaten';

data = dlmread('AtmosphericData_July_fs10Hz_Kurz.txt','\t');
data_name = 'athmosphärische Daten';

%overall mean/std 
mean = nanmean(data);
std = nanstd(data);
mag_fluc = std^2;
degree_turb = mag_fluc/mean^2;
u_dash_tot = data-mean;
u_dash_tot_norm = u_dash_tot/std;

lines = size(data,1);

%time men/std
t = 6000;
for i = 1 : 1:lines/t
    mean_t1(i) =nanmean (data(i*t-(t-1):i*t)); 
    std_t1(i) =nanstd (data(i*t-(t-1):i*t)); 
    mag_fluc_t1(i) = std_t1(i)^2;
    u_dash_t1Mean(i*t-(t-1):i*t) = data(i*t-(t-1):i*t)- mean_t1(i);
    u_dash_t1Mean_norm =  u_dash_t1Mean(i*t-(t-1):i*t)./std_t1(i);
    
end

degree_turb_t1 = mag_fluc_t1./mean_t1.^2;
clear t;

t = 600;
for i = 1 : 1:lines/t
    mean_t2(i) =nanmean (data(i*t-(t-1):i*t)); 
    std_t2(i) =nanstd (data(i*t-(t-1):i*t));
    mag_fluc_t2(i) = std_t2(i)^2;    
    u_dash_t2Mean(i*t-(t-1):i*t) = data(i*t-(t-1):i*t)- mean_t2(i);  
    u_dash_t2Mean_norm =  u_dash_t2Mean(i*t-(t-1):i*t)./std_t2(i);
    
end
degree_turb_t2 = mag_fluc_t2./mean_t2.^2;
clear t;

t = 60;
for i = 1 : 1:lines/t
    mean_t3(i) =nanmean (data(i*t-(t-1):i*t)); 
    std_t3(i) =nanstd (data(i*t-(t-1):i*t));
    mag_fluc_t3(i) = std_t3(i)^2;
    u_dash_t3Mean(i*t-(t-1):i*t) = data(i*t-(t-1):i*t)- mean_t3(i);
    u_dash_t3Mean_norm =  u_dash_t3Mean(i*t-(t-1):i*t)./std_t3(i);
    
end
degree_turb_t3 = mag_fluc_t3./mean_t3.^2;
clear t;


%table: Mean Value
lines_name_tab = {'average [m/s]','standard deviation [m/s]'};
mean_6000 = [nanmean(mean_t1);nanstd(mean_t1)];
mean_600 = [nanmean(mean_t2);nanstd(mean_t2)];
mean_60 = [nanmean(mean_t3);nanstd(mean_t3)];
mean_comp = [mean;0];

Tabele_mean= table(mean_6000,mean_600,mean_60,mean_comp,'RowNames',lines_name_tab);

%table: magnitude of fluctuation
lines_name_tab = {'average [m^2/s^2]','standard deviation [m^2/s^2]'};
mag_fluc_6000 = [nanmean(mag_fluc_t1);nanstd(mag_fluc_t1)];
mag_fluc_600 = [nanmean(mag_fluc_t2);nanstd(mag_fluc_t2)];
mag_fluc_60 = [nanmean(mag_fluc_t3);nanstd(mag_fluc_t3)];
mag_fluc_comp = [mag_fluc;0];

Tabele_mag_fluc= table(mag_fluc_6000,mag_fluc_600,mag_fluc_60,mag_fluc_comp,'RowNames',lines_name_tab);


%table: Degree of turbulence
lines_name_tab = {'average [1]','standard deviation [1]'};
degree_of_turb_6000 = [nanmean(degree_turb_t1);nanstd(degree_turb_t1)];
degree_of_turb_600 = [nanmean(degree_turb_t2);nanstd(degree_turb_t2)];
degree_of_turb_60 = [nanmean(degree_turb_t3);nanstd(degree_turb_t3)];
degree_of_turb_comp = [degree_turb;0];

Tabele_degree_turb= table(degree_of_turb_6000,degree_of_turb_600,degree_of_turb_60,degree_of_turb_comp,'RowNames',lines_name_tab);

%PDF,u
    
figure
hold on
h = histogram(data,40,'Normalization','probability');
title(['Relative Häufigkeit von u, ',data_name], 'FontSize', 15)
xlabel('u [m/s]', 'FontSize', 15)
ylabel('Relative Häufigkeit [1]', 'FontSize', 15)
hold off
saveas (gcf,['PDF of u, ',data_name,'.png']);



figure
hold on
h = histogram(u_dash_tot,40,'Normalization','probability','FaceColor','r')
h = histogram(u_dash_t1Mean,40,'Normalization','probability','FaceColor','b')
h = histogram(u_dash_t2Mean,40,'Normalization','probability','FaceColor','g')
h = histogram(u_dash_t3Mean,40,'Normalization','probability','FaceColor','y')
title(['Relative Häufigkeit von u\prime, ',data_name], 'FontSize', 15)
xlabel('u\prime [m/s]', 'FontSize', 15)
ylabel('Relative Häufigkeit [1]', 'FontSize', 15)
legend('Gesamter Datensatz','6000 Daten Intervall','600 Daten Intervall','60 Daten Intervall')
hold off
saveas (gcf,['PDF of u_dash, ',data_name,'.png']);


figure
hold on
h = histogram(u_dash_tot_norm,80,'Normalization','probability','FaceColor','r')
h = histogram(u_dash_t1Mean_norm,80,'Normalization','probability','FaceColor','b')
h = histogram(u_dash_t2Mean_norm,80,'Normalization','probability','FaceColor','g')
h = histogram(u_dash_t3Mean_norm,80,'Normalization','probability','FaceColor','y')
title(['Relative Häufigkeit von u\prime/\sigma, ',data_name], 'FontSize', 15)
xlabel('u\prime/\sigma [1]', 'FontSize', 15)
ylabel('Relative Häufigkeit [1]', 'FontSize', 15)
legend('Gesamter Datensatz','6000 Daten Intervall','600 Daten Intervall','60 Daten Intervall')
hold off
saveas (gcf,['PDF of u_dash_norm, ',data_name,'.png'])