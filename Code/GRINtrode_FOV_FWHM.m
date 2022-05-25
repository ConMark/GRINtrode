function GRINtrode_FOV_FWHM
 
data = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\fov_max.csv');
figNo=0;
 
pos = transpose(data(:,1));
val = transpose(data(:,2));
f = fit(pos.',val.','gauss2')
 
figNo=figNo+1;
figure(figNo)
plot(pos,val,'.-','LineWidth',4,'MarkerSize',32)
hold on
fitplot=plot(f)
set(fitplot,'lineWidth',4);
xlim([0 max(pos)])
ylim([750 1800])
legend('Data','Fitted Gaussian','FontSize',44)
title('FOV field profile for imaging system A', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 60)
ylabel('Intensity', 'FontSize', 60)
ax = gca;
ax.FontSize = 54; 
 
FWHM = 2*sqrt(log(2))*f.c1
 
 
 


