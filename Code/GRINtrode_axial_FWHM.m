function GRINtrode_axial_FWHM
 
data1 = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\ax1.csv');
data2 = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\ax2.csv');
data3 = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\ax3.csv');
data4 = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\ax4.csv');
data5 = readmatrix('C:\Users\Connor\Documents\Work\Paper\GT9mmGT14mmNewResFOV\take2\ax5.csv');
 
pos1 = transpose(data1(:,1));
val1 = transpose(data1(:,2));
pos2 = transpose(data2(:,1));
val2 = transpose(data2(:,2));
pos3 = transpose(data3(:,1));
val3 = transpose(data3(:,2));
pos4 = transpose(data4(:,1));
val4 = transpose(data4(:,2));
pos5 = transpose(data5(:,1));
val5 = transpose(data5(:,2));
 
f1 = fit(pos1.',val1.','gauss2')
f2 = fit(pos2.',val2.','gauss2')
f3 = fit(pos3.',val3.','gauss2')
f4 = fit(pos4.',val4.','gauss2')
f5 = fit(pos5.',val5.','gauss2')
 
FWHM1 = 2*sqrt(log(2))*f1.c1
FWHM2 = 2*sqrt(log(2))*f2.c1
FWHM3 = 2*sqrt(log(2))*f3.c1
FWHM4 = 2*sqrt(log(2))*f4.c1
FWHM5 = 2*sqrt(log(2))*f5.c1
 
avgFWHM=mean([FWHM1,FWHM2,FWHM3,FWHM4,FWHM5])
 
figNo=0;
 
figNo=figNo+1;
figure(figNo)
plot(pos1,val1,'.-','LineWidth',2,'MarkerSize',16)
hold on
fitplot=plot(f1);
set(fitplot,'lineWidth',2);
xlim([0 max(pos1)])
% ylim([750 1800])
legend('Data','Fitted Gaussian','FontSize',36)
title('Example 0.5 micron bead line profile for imaging system B', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 44)
ylabel('Intensity', 'FontSize', 44)
ax = gca;
ax.FontSize = 36; 
 
figNo=figNo+1;
figure(figNo)
plot(pos2,val2,'.-','LineWidth',2,'MarkerSize',16)
hold on
fitplot=plot(f2);
set(fitplot,'lineWidth',2);
xlim([0 max(pos2)])
% ylim([750 1800])
legend('Data','Fitted Gaussian','FontSize',36)
title('Example 0.5 micron bead line profile for imaging system B', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 44)
ylabel('Intensity', 'FontSize', 44)
ax = gca;
ax.FontSize = 36; 
 
figNo=figNo+1;
figure(figNo)
plot(pos3,val3,'.-','LineWidth',2,'MarkerSize',16)
hold on
fitplot=plot(f3);
set(fitplot,'lineWidth',2);
xlim([0 max(pos3)])
% ylim([750 1800])
legend('Data','Fitted Gaussian','FontSize',36)
title('Example 0.5 micron bead axial profile for imaging system B', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 44)
ylabel('Intensity', 'FontSize', 44)
ax = gca;
ax.FontSize = 36; 
 
figNo=figNo+1;
figure(figNo)
plot(pos4,val4,'.-','LineWidth',2,'MarkerSize',16)
hold on
fitplot=plot(f4);
set(fitplot,'lineWidth',2);
xlim([0 max(pos4)])
% ylim([750 1800])
legend('Data','Fitted Gaussian','FontSize',36)
title('Example 0.5 micron bead axial profile for imaging system B', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 44)
ylabel('Intensity', 'FontSize', 44)
ax = gca;
ax.FontSize = 36; 
 
figNo=figNo+1;
figure(figNo)
plot(pos5,val5,'.-','LineWidth',4,'MarkerSize',32)
hold on
fitplot=plot(f5);
set(fitplot,'lineWidth',4);
xlim([0 max(pos5)])
ylim([200 450])
legend('Data','Fitted Gaussian','FontSize',44)
title('Example 0.5 micron bead axial profile for imaging system B', 'FontSize', 72)
xlabel('Position (microns)', 'FontSize', 60)
ylabel('Intensity', 'FontSize', 60)
ax = gca;
ax.FontSize = 54; 

