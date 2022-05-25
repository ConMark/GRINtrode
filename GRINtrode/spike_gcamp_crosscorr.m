 function spike_gcamp_crosscorr

%run read_Intan_RHD2000_file.m and drag in jt times before running

T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii

optical_data = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\DFTemporal.csv');
optical_data = optical_data/100;
spatial = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\Spatial.csv'); %cell CoM locations from CaImAn
avgimg = imread('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\15765_9_27_TL2_MAX_denoised_scalebar.png'); %max projection of timelapse
framerate = 3.848522167; %framerate in frames/second
decay_constant = 1.56; %decay constant from caiman
save_results = 1; %save_results = 1 will save all matrices and figures

if save_results == 1
    t = datetime('now');
    t = datestr(t);
    t = strrep(t,':','_');
    t = strrep(t,' ','_');
    t = strrep(t,'-','_');
    folder = ['C:\Users\Connor\Desktop\157565_9_27_2019_TL2_spike_gcamp_crosscorr_Diegos_waveclus_', t, '\'];
    mkdir(folder);
end

intan_freq = frequency_parameters.amplifier_sample_rate %intan sampling rate in Hz
trig=find(board_dig_in_data(6,:),1,'first') %find trigger input from imaging start
imaging_start_time = trig*(1/intan_freq) %when imaging started after starting intan recording
imaging_length_seconds = size(optical_data,2)/framerate %length of timelapse
intan_points = int64(imaging_length_seconds*20000) %number of intan data points

imaging_time = zeros(1,size(optical_data,2));
for i=1:size(optical_data,2)
    imaging_time(1,i) = i/framerate;
end

figNo=0;
offsets = offset_for_chan;
timestamps = all_timestamp_per_file;
spikes = cluster_class_per_file;
spikes = spikes+1;

%tetrode 1
tetrode_1_clusters = spikes(1:offsets(2));
tetrode_1_times = timestamps(1:offsets(2));
tetrode_1_indexes = tetrode_1_times*20000;
max_timestamps(1) = max(tetrode_1_times);

%tetrode 2
tetrode_2_clusters = spikes(offsets(2)+1:offsets(3));
tetrode_2_times = timestamps(offsets(2)+1:offsets(3));
tetrode_2_indexes = tetrode_2_times*20000;
max_timestamps(2) = max(tetrode_2_times);

%tetrode 3
tetrode_3_clusters = spikes(offsets(3)+1:offsets(4));
tetrode_3_times = timestamps(offsets(3)+1:offsets(4));
tetrode_3_indexes = tetrode_3_times*20000;
max_timestamps(3) = max(tetrode_3_times);

%tetrode 4
tetrode_4_clusters = spikes(offsets(4)+1:size(spikes,2));
tetrode_4_times = timestamps(offsets(4)+1:size(spikes,2));
tetrode_4_indexes = tetrode_4_times*20000;
max_timestamps(4) = max(tetrode_4_times);

%%generate time 
time = zeros(1,intan_points);
for i = 0:intan_points-1
    time(i+1) = i/20000;
end

%%tetrode 1 spikes
num_cells1 = max(tetrode_1_clusters(1,:))-1
tetrode_1_spikes = zeros(num_cells1, size(amplifier_data,2));
tetrode_1_indexes = int64(tetrode_1_indexes);
for i = 1:size(tetrode_1_clusters,2)
    tetrode_1_spikes(tetrode_1_clusters(i),tetrode_1_indexes(i)) = 1;
end

tetrode_1_spikes = tetrode_1_spikes(:,[trig:trig+intan_points-1]);

figNo = figNo+1;
figure(figNo)
for i = 1:num_cells1
    plot(time, tetrode_1_spikes(i+1,:)+2*(i-1))
    hold on
end
title('Tetrode 1 Clusters')

%%tetrode 2 spikes
num_cells2 = max(tetrode_2_clusters(1,:))-1
tetrode_2_spikes = zeros(num_cells2, size(amplifier_data,2));
tetrode_2_indexes = int64(tetrode_2_indexes);
for i = 1:size(tetrode_2_clusters,2)
    tetrode_2_spikes(tetrode_2_clusters(i),tetrode_2_indexes(i)) = 1;
end

tetrode_2_spikes = tetrode_2_spikes(:,[trig:trig+intan_points-1]);

figNo = figNo+1;
figure(figNo)
for i = 1:num_cells2
    plot(time, tetrode_2_spikes(i+1,:)+2*(i-1))
    hold on
end
title('Tetrode 2 Clusters')

%%tetrode 3 spikes
num_cells3 = max(tetrode_3_clusters(1,:))-1
tetrode_3_spikes = zeros(num_cells3, size(amplifier_data,2));
tetrode_3_indexes = int64(tetrode_3_indexes);
for i = 1:size(tetrode_3_clusters,2)
    tetrode_3_spikes(tetrode_3_clusters(i),tetrode_3_indexes(i)) = 1;
end

tetrode_3_spikes = tetrode_3_spikes(:,[trig:trig+intan_points-1]);

figNo = figNo+1;
figure(figNo)
for i = 1:num_cells3
    plot(time, tetrode_3_spikes(i+1,:)+2*(i-1),'LineWidth',2)
    hold on
end
xlim([0 max(time)])
% ylim([0 max(optical_data(max(cells),:)+100*(j-2))+10])
ylim([0 3.5])

ax = gca;
set(gca,'YTick',[])
ax.FontSize = 36; 
title('Single Unit Spikes', 'FontSize', 56)
xlabel('Time (s)', 'FontSize', 44)


%%tetrode 4 spikes
num_cells4 = max(tetrode_4_clusters(1,:))-1
tetrode_4_spikes = zeros(num_cells4, size(amplifier_data,2));
tetrode_4_indexes = int64(tetrode_4_indexes);
for i = 1:size(tetrode_4_clusters,2)
    tetrode_4_spikes(tetrode_4_clusters(i),tetrode_4_indexes(i)) = 1;
end

tetrode_4_spikes = tetrode_4_spikes(:,[trig:trig+intan_points-1]);

figNo = figNo+1;
figure(figNo)
for i = 1:num_cells4
    plot(time, tetrode_4_spikes(i+1,:)+2*(i-2))
    hold on
end
title('Tetrode 4 Clusters')

%exponential decay function
for i =1:20000*10;
    g(i) = exp(-time(i)/decay_constant);
end

figNo = figNo+1;
figure(figNo)
plot(time(1:20000*10),g, 'LineWidth', 3)
title('Exponential Decay Function','FontSize', 52)

figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 1 raw (top) vs convolved (bottom)')
for i = 2:num_cells1+1
    tetrode_1_convolved1(i-1,:) = conv(g,tetrode_1_spikes(i,:));
    tetrode_1_convolved(i-1,:) = tetrode_1_convolved1(i-1, [1:intan_points]);
    subplot(2,num_cells1,i-1)
    plot(time,tetrode_1_spikes(i,:))
    subplot(2,num_cells1,i-1+num_cells1)
    plot(time, tetrode_1_convolved(i-1,:))
end



figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 2 raw (top) vs convolved (bottom)')
for i = 2:num_cells2+1
    tetrode_2_convolved1(i-1,:) = conv(g,tetrode_2_spikes(i,:));
    tetrode_2_convolved(i-1,:) = tetrode_2_convolved1(i-1, [1:intan_points]);
    subplot(2,num_cells2,i-1)
    plot(time,tetrode_2_spikes(i,:))
    subplot(2,num_cells2,i-1+num_cells2)
    plot(time, tetrode_2_convolved(i-1,:))
end


figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 3 raw (top) vs convolved (bottom)')
for i = 2:num_cells3+1
    tetrode_3_convolved1(i-1,:) = conv(g,tetrode_3_spikes(i,:));
    tetrode_3_convolved(i-1,:) = tetrode_3_convolved1(i-1, [1:intan_points]);
    subplot(2,num_cells3,i-1)
    plot(time,tetrode_3_spikes(i,:))
    subplot(2,num_cells3,i-1+num_cells3)
    plot(time, tetrode_3_convolved(i-1,:))
end


figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 4 raw (top) vs convolved (bottom)')
for i = 2:num_cells4+1
    tetrode_4_convolved1(i-1,:) = conv(g,tetrode_4_spikes(i,:));
    tetrode_4_convolved(i-1,:) = tetrode_4_convolved1(i-1, [1:intan_points]);
    subplot(2,num_cells4,i-1)
    plot(time,tetrode_4_spikes(i,:))
    subplot(2,num_cells4,i-1+num_cells4)
    plot(time, tetrode_4_convolved(i-1,:))
end


dec_n=100;
dectim=downsample(time,dec_n);

for i=1:size(optical_data,1)
    DFF_interp(i,:)=interp1(imaging_time, optical_data(i,:), dectim, 'linear','extrap');
end

for i = 1:size(tetrode_1_convolved,2)
    convtime(i) = i/20000;
end

for i = 1:num_cells1
    decspikes1(i,:) = interp1(convtime,tetrode_1_convolved(i,:),dectim, 'linear', 'extrap');
end

for i = 1:num_cells2
    decspikes2(i,:) = interp1(convtime,tetrode_2_convolved(i,:),dectim, 'linear', 'extrap');
end

for i = 1:num_cells3
    decspikes3(i,:) = interp1(convtime,tetrode_3_convolved(i,:),dectim, 'linear', 'extrap');
end

for i = 1:num_cells4
    decspikes4(i,:) = interp1(convtime,tetrode_4_convolved(i,:),dectim, 'linear', 'extrap');
end

size(dectim)
size(decspikes1)
size(DFF_interp)

figNo=figNo+1
figure(figNo)
sgtitle('Tetrode 1 convolved (top) vs convolved decimated (bottom)')
for i = 1:num_cells1
    subplot(2,num_cells1,i)
    plot(convtime,tetrode_1_convolved(i,:))
    subplot(2,num_cells1,i+num_cells1)
    plot(dectim,decspikes1(i,:))
end


figNo=figNo+1
figure(figNo)
sgtitle('Tetrode 2 convolved (top) vs convolved decimated (bottom)')
for i = 1:num_cells2
    subplot(2,num_cells2,i)
    plot(convtime,tetrode_2_convolved(i,:))
    subplot(2,num_cells2,i+num_cells2)
    plot(dectim,decspikes2(i,:))
end


figNo=figNo+1
figure(figNo)
sgtitle('Tetrode 3 convolved (top) vs convolved decimated (bottom)')
for i = 1:num_cells3
    subplot(2,num_cells3,i)
    plot(convtime,tetrode_3_convolved(i,:))
    subplot(2,num_cells3,i+num_cells3)
    plot(dectim,decspikes3(i,:))
end


figNo=figNo+1
figure(figNo)
sgtitle('Tetrode 4 convolved (top) vs convolved decimated (bottom)')
for i = 1:num_cells4
    subplot(2,num_cells4,i)
    plot(convtime,tetrode_4_convolved(i,:))
    subplot(2,num_cells4,i+num_cells4)
    plot(dectim,decspikes4(i,:))
end


for i = 1:size(optical_data,1)
    for j =1:num_cells1
    [corr,P] = corrcoef(decspikes1(j,:),DFF_interp(i,:));
    corrvals1(i,j) = corr(1,2);
    P1(i,j) = P(1,2);
    end
end

for i = 1:num_cells1
    [HighestCorr1(i), CellNumber1(i)] = max(corrvals1(:,i));
    MaxP1(i) = P1(CellNumber1(i), i);
end

HighestCorr1
CellNumber1


for i = 1:size(optical_data,1)
    for j =1:num_cells2
    [corr,P] = corrcoef(decspikes2(j,:),DFF_interp(i,:));
    corrvals2(i,j) = corr(1,2);
    P2(i,j) = P(1,2);
    end
end

for i = 1:num_cells2
    [HighestCorr2(i), CellNumber2(i)] = max(corrvals2(:,i));
    MaxP2(i) = P2(CellNumber2(i), i);
end
[HighestCorr2Sorted,sortIdx] = sort(corrvals2(:,2),'descend');
HighestCorr2Sorted
cells=[1:1:size(optical_data,2)];
CellNumber2sorted = cells(sortIdx)

HighestCorr2
CellNumber2

for i = 1:size(optical_data,1)
    for j =1:num_cells3
    [corr,P] = corrcoef(decspikes3(j,:),DFF_interp(i,:));
    corrvals3(i,j) = corr(1,2);
    P3(i,j) = P(1,2);
    end
end

for i = 1:num_cells3
    [HighestCorr3(i), CellNumber3(i)] = max(corrvals3(:,i));
    MaxP3(i) = P3(CellNumber3(i), i);
end

HighestCorr3
CellNumber3


for i = 1:size(optical_data,1)
    for j =1:num_cells4
    [corr,P] = corrcoef(decspikes4(j,:),DFF_interp(i,:));
    corrvals4(i,j) = corr(1,2);
    P4(i,j) = P(1,2);
    end
end

for i = 1:num_cells4
    [HighestCorr4(i), CellNumber4(i)] = max(corrvals4(:,i));
    MaxP4(i) = P4(CellNumber4(i), i);
end

HighestCorr4
CellNumber4

figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 1 convolved & decimated clusters (top) vs interpolated highest correlation DFF (bottom)')
for i = 1:num_cells1
    subplot(2,num_cells1,i)
    plot(dectim,decspikes1(i,:))
    title("Tetrode 1, cluster number " + i)
    subplot(2,num_cells1,i+num_cells1)
    plot(dectim,DFF_interp(CellNumber1(i),:))
    title("Cell number = " + CellNumber1(i) + ", rho = " + HighestCorr1(i))
end


figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 2 convolved & decimated clusters (top) vs interpolated highest correlation DFF (bottom)')
for i = 1:num_cells2
    subplot(2,num_cells2,i)
    plot(dectim,decspikes2(i,:))
    title("Tetrode 2, cluster number " + i)
    subplot(2,num_cells2,i+num_cells2)
    plot(dectim,DFF_interp(CellNumber2(i),:))
    title("Cell number = " + CellNumber2(i) + ", rho = " + HighestCorr2(i))
end


figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 3 convolved & decimated clusters (top) vs interpolated highest correlation DFF (bottom)')
for i = 1:num_cells3
    subplot(2,num_cells3,i)
    plot(dectim,decspikes3(i,:))
    title("Tetrode 3, cluster number " + i)
    subplot(2,num_cells3,i+num_cells3)
    plot(dectim,DFF_interp(CellNumber3(i),:))
    title("Cell number = " + CellNumber3(i) + ", rho = " + HighestCorr3(i))
end


figNo = figNo+1;
figure(figNo)
sgtitle('Tetrode 4 convolved & decimated clusters (top) vs interpolated highest correlation DFF (bottom)')
for i = 1:num_cells4
    subplot(2,num_cells4,i)
    plot(dectim,decspikes4(i,:))
    title("Tetrode 4, cluster number " + i)
    subplot(2,num_cells4,i+num_cells4)
    plot(dectim,DFF_interp(CellNumber4(i),:))
    title("Cell number = " + CellNumber4(i) + ", rho = " + HighestCorr4(i))
end

LabeledImage = insertMarker(avgimg,[spatial(CellNumber1(:),2) spatial(CellNumber1(:),1)],'color','red');
figNo=figNo+1;
figure(figNo)
imshow(LabeledImage)
title('tetrode 1 cells')

LabeledImage = insertMarker(avgimg,[spatial(CellNumber2(:),2) spatial(CellNumber2(:),1)],'color','green');
figNo=figNo+1;
figure(figNo)
imshow(LabeledImage)
title('tetrode 2 cells')

LabeledImage = insertMarker(avgimg,[spatial(CellNumber3(:),2) spatial(CellNumber3(:),1)],'color','blue');
figNo=figNo+1;
figure(figNo)
imshow(LabeledImage)
title('tetrode 3 cells')

LabeledImage = insertMarker(avgimg,[spatial(CellNumber4(:),2) spatial(CellNumber4(:),1)],'color','magenta');
figNo=figNo+1;
figure(figNo)
imshow(LabeledImage)
title('tetrode 4 cells')

figNo = figNo+1;
figure(figNo)
plot(dectim,decspikes2(2,:)*(1/3.57),'r','LineWidth',2)
hold on
plot(dectim,DFF_interp(CellNumber2(2),:)-0.7729,'g','LineWidth',2)
% hold on
% plot(dectim,DFF_interp(59,:),'b','LineWidth',2)
% hold on
% plot(dectim,DFF_interp(77,:),'b','LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('Cross-Correlation of GCaMP and Convolved Spike Data, \rho = 0.48', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'Spike Data', 'GCaMP Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo = figNo+1;
figure(figNo)
plot(dectim,decspikes3(2,:),'r','LineWidth',2)
hold on
plot(dectim,DFF_interp(CellNumber3(2),:)-5.4,'b','LineWidth',2)
% hold on
% plot(dectim,DFF_interp(59,:),'b','LineWidth',2)
% hold on
% plot(dectim,DFF_interp(77,:),'b','LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('Cross-Correlation of GCaMP and Convolved Spike Data, \rho = ', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'Spike Data', 'GCaMP Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

% mdl=fitglm(dFFtraces',sgf_this_decimated_LFP_logPtraces,'linear');
% LFPlog_pred = predict(mdl,dFFtraces');
% 
% order=3;
% framelen=31;
% sgf_this_decimated_LFP_logPtraces = sgolayfilt(this_decimated_LFP_logPtraces,order,framelen);
% mdl=fitglm(dFFtraces',sgf_this_decimated_LFP_logPtraces,'linear');
% LFPlog_pred = predict(mdl,dFFtraces');

mdl11=fitglm(DFF_interp',decspikes1(1,:),'linear');
LFPlog_pred11 = predict(mdl11,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred11,decspikes1(1,:));
corrvals11 = corr(1,2)

mdl12=fitglm(DFF_interp',decspikes1(2,:),'linear');
LFPlog_pred12 = predict(mdl12,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred12,decspikes1(2,:));
corrvals12 = corr(1,2)

mdl21=fitglm(DFF_interp',decspikes2(1,:),'linear');
LFPlog_pred21 = predict(mdl21,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred21,decspikes2(1,:));
corrvals21 = corr(1,2)

mdl22=fitglm(DFF_interp',decspikes2(2,:),'linear');
LFPlog_pred22 = predict(mdl22,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred22,decspikes2(2,:));
corrvals22 = corr(1,2)

mdl31=fitglm(DFF_interp',decspikes3(1,:),'linear');
LFPlog_pred31 = predict(mdl31,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred31,decspikes3(1,:));
corrvals31 = corr(1,2)

mdl32=fitglm(DFF_interp',decspikes3(2,:),'linear');
LFPlog_pred32 = predict(mdl32,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred32,decspikes3(2,:));
corrvals32 = corr(1,2)

mdl41=fitglm(DFF_interp',decspikes4(1,:),'linear');
LFPlog_pred41 = predict(mdl41,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred41,decspikes4(1,:));
corrvals41 = corr(1,2)

mdl42=fitglm(DFF_interp',decspikes4(2,:),'linear');
LFPlog_pred42 = predict(mdl42,DFF_interp');

[corr,P] = corrcoef(LFPlog_pred42,decspikes4(2,:));
corrvals42 = corr(1,2)

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred11,'LineWidth',2)
hold on
plot(dectim,decspikes1(1,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 1,1 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred12,'LineWidth',2)
hold on
plot(dectim,decspikes1(2,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 1,2 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred21,'LineWidth',2)
hold on
plot(dectim,decspikes2(1,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 2,1 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred22,'LineWidth',2)
hold on
plot(dectim,decspikes2(2,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 2,2 Results of GCaMP and Convolved Spike Data, \rho = ', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred31,'LineWidth',2)
hold on
plot(dectim,decspikes3(1,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 3,1 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred32,'LineWidth',2)
hold on
plot(dectim,decspikes3(2,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 3,2 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred41,'LineWidth',2)
hold on
plot(dectim,decspikes4(1,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 4,1 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])

figNo=figNo+1;
figure(figNo)
plot(dectim,LFPlog_pred42,'LineWidth',2)
hold on
plot(dectim,decspikes4(2,:),'LineWidth',2)
xlabel('Time (sec)', 'FontSize', 38)
ylabel('\DeltaF/F (GCaMP Data) / Arbitrary (Spike Data)', 'FontSize', 38);
title('GLM 4,2 Results of GCaMP and Convolved Spike Data', 'FontSize', 52)
ax = gca;
ax.FontSize = 28;
legend({'GCaMP Data', 'Spike Data'}, 'FontSize', 28)
xlim([0 max(dectim)])


if save_results == 1
    %save all figures
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = num2str(get(FigHandle, 'Number'));
      set(0, 'CurrentFigure', FigHandle);
      savefig(fullfile(folder, [FigName '.fig']));
    end
    %write matrices
    writematrix(HighestCorr1, strcat(folder,'HighestCorr1.txt'))
    writematrix(CellNumber1, strcat(folder,'CellNumber1.txt'))
    writematrix(HighestCorr2, strcat(folder,'HighestCorr2.txt'))
    writematrix(CellNumber2, strcat(folder,'CellNumber2.txt'))
    writematrix(HighestCorr3, strcat(folder,'HighestCorr3.txt'))
    writematrix(CellNumber3, strcat(folder,'CellNumber3.txt'))
    writematrix(HighestCorr4, strcat(folder,'HighestCorr4.txt'))
    writematrix(CellNumber4, strcat(folder,'CellNumber4.txt'))
end