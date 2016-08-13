clear all
close all

%Import data
data_lab = dlmread('jfm_data_block1.txt','\t');


%overall mean/std 
mean_lab = nanmean(data_lab);
std_lab = nanstd(data_lab);
degree_turb_lab = std_lab.^2/mean_lab.^2;

lines = size(data_lab,1);

%time men/std
t = 6000;
for i = 1 : 1:lines/t
    mean_t1_lab(i) =nanmean (data_lab(i*t-(t-1):i*t)); 
    std_t1_lab(i) =nanstd (data_lab(i*t-(t-1):i*t)); 
    
    u_dash_t1Mean(i*t-(t-1):i*t) = data_lab(i*t-(t-1):i*t)- mean_t1_lab(i);
    
    u_dash_t1Mean_norm =  u_dash_t1Mean(i*t-(t-1):i*t)./std_t1_lab(i);
    
end
degree_turb_t1_lab = std_t1_lab.^2./mean_t1_lab.^2;
clear t;

t = 600;
for i = 1 : 1:lines/t
    mean_t2_lab(i) =nanmean (data_lab(i*t-(t-1):i*t)); 
    std_t2_lab(i) =nanstd (data_lab(i*t-(t-1):i*t)); 
    
    u_dash_t2Mean(i*t-(t-1):i*t) = data_lab(i*t-(t-1):i*t)- mean_t2_lab(i);
    
     u_dash_t2Mean_norm =  u_dash_t2Mean(i*t-(t-1):i*t)./std_t2_lab(i);
    
end
degree_turb_t2_lab = std_t2_lab.^2./mean_t2_lab.^2;
clear t;

t = 60;
for i = 1 : 1:lines/t
    mean_t3_lab(i) =nanmean (data_lab(i*t-(t-1):i*t)); 
    std_t3_lab(i) =nanstd (data_lab(i*t-(t-1):i*t)); 
    
    u_dash_t3Mean(i*t-(t-1):i*t) = data_lab(i*t-(t-1):i*t)- mean_t3_lab(i);
    
     u_dash_t3Mean_norm =  u_dash_t3Mean(i*t-(t-1):i*t)./std_t3_lab(i);
    
end
degree_turb_t3_lab = std_t3_lab.^2./mean_t3_lab.^2;
clear t;

%table: Mean Value
lines_name_tab = {'min','max'};
std_6000 = [min(std_t1_lab);max(std_t1_lab)];
std_600 = [min(std_t2_lab);max(std_t2_lab)];
std_60 = [min(std_t3_lab);max(std_t3_lab)];
std_comp = [std_lab;std_lab];

Tabele_std= table(std_6000,std_600,std_60,std_comp,'RowNames',lines_name_tab);


%table: standart deviation
lines_name_tab = {'min','max'};
mean_6000 = [min(mean_t1_lab);max(mean_t1_lab)];
mean_600 = [min(mean_t2_lab);max(mean_t2_lab)];
mean_60 = [min(mean_t3_lab);max(mean_t3_lab)];
mean_comp = [mean_lab;mean_lab];

Tabele_mean= table(mean_6000,mean_600,mean_60,mean_comp,'RowNames',lines_name_tab);

%table: Degree of turbulence
lines_name_tab = {'min','max'};
degree_of_turb_6000 = [min(degree_turb_t1_lab);max(degree_turb_t1_lab)];
degree_of_turb_600 = [min(degree_turb_t2_lab);max(degree_turb_t2_lab)];
degree_of_turb_60 = [min(degree_turb_t3_lab);max(degree_turb_t3_lab)];
degree_of_turb_comp = [degree_turb_lab;degree_turb_lab];

Tabele_degree_turb= table(degree_of_turb_6000,degree_of_turb_600,degree_of_turb_60,degree_of_turb_comp,'RowNames',lines_name_tab);

%PDF,u
%pts = 0:0.1:5;
[f_u_tot,xi_utot] =  ksdensity (data_lab);     
figure
plot (xi_utot,f_u_tot)
title('PDF, data from laboratoty')
xlabel('wind speed [m/s]')
ylabel('probability [1]')
clear xi;
clear f;

%PDF,u_dash

u_dash_totalMean = data_lab-mean_lab;

[f_dash_total,xi_u_dash_total] =  ksdensity (u_dash_totalMean); 
[f_dash_t1,xi_u_dash_t1] =  ksdensity (u_dash_t1Mean); 
[f_dash_t2,xi_u_dash_t2] =  ksdensity (u_dash_t2Mean);
[f_dash_t3,xi_u_dash_t3] =  ksdensity (u_dash_t3Mean);

figure
hold on
plot (xi_u_dash_t1,f_dash_t1)
plot (xi_u_dash_t2,f_dash_t2)
plot (xi_u_dash_t3,f_dash_t3)
plot(xi_u_dash_total,f_dash_total)
hold off

title('PDF fluctuational part of velocity, data from laboratoty')
xlabel('wind speed fluctuation [m/s]')
ylabel('probability [1]')
legend('averaging intervall = 6000','averaging intervall = 600','averaging intervall = 60','total average')

%PDF,u_dash, normalisiert y its standartdeviation

u_dash_totalMean_norm = (data_lab-mean_lab)./std_lab;

[f_dash_total_norm,xi_u_dash_total_norm] =  ksdensity (u_dash_totalMean); 
[f_dash_t1_norm,xi_u_dash_t1_norm] =  ksdensity (u_dash_t1Mean_norm); 
[f_dash_t2_norm,xi_u_dash_t2_norm] =  ksdensity (u_dash_t2Mean_norm);
[f_dash_t3_norm,xi_u_dash_t3_norm] =  ksdensity (u_dash_t3Mean_norm);

figure
hold on
plot (xi_u_dash_t1_norm,f_dash_t1_norm)
plot (xi_u_dash_t2_norm,f_dash_t2_norm)
plot (xi_u_dash_t3_norm,f_dash_t3_norm)
plot(xi_u_dash_total_norm,f_dash_total_norm)
hold off

title('PDF fluctuational part of velocity, data from laboratoty')
xlabel('wind speed fluctuation, normalized by standart deviation [1]')
ylabel('probability [1]')
legend('averaging intervall = 6000','averaging intervall = 600','averaging intervall = 60','total average')



