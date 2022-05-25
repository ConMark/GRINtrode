%drgcmCrossPCArhddFFv2
%clear all
close all

rng('shuffle')

use_raw=1; %If this is 1 the program uses the raw trace, otherwise it uses the inferred trace

framerate = 3.848522167;
no_files=1;
no_electrodes=16;

T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii

optical_data = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\DFTemporal.csv'); %DF/F traces from CaImAn
%optical_data = optical_data/100;
spatial = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\Spatial.csv'); %cell CoM locations from CaImAn
avgimg = imread('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\15765_9_27_TL2_MAX_denoised_scalebar.png'); %max projection of timelapse
framerate = 3.848522167; %framerate in frames/second
save_results = 1;
maxfreq=100;
no_electrodes=16;

size(optical_data)

if save_results == 1
    t = datetime('now')
    t = datestr(t)
    t = strrep(t,':','_')
    t = strrep(t,' ','_')
    t = strrep(t,'-','_')
    folder = ['C:\Users\Connor\Desktop\DanielAnimal9_TL12_2_25_2022_bw_power_crosscor_', t, '\']
    mkdir(folder);
end

intan_freq = frequency_parameters.amplifier_sample_rate %intan sampling rate in Hz
trig=find(board_dig_in_data(6,:),1,'first'); %find trigger input from imaging start
imaging_start_time = trig*(1/intan_freq) %when imaging started after starting intan recording
imaging_length_seconds = size(optical_data,2)/framerate %length of timelapse
intan_points = int64(imaging_length_seconds*20000) %number of intan data points

%create tetrode data array
tetrode_data = zeros(no_electrodes,intan_points);
for i = 1:no_electrodes
    tetrode_data(i,:) = amplifier_data(i+8,[trig:trig+intan_points-1]); %for data with electrodes on channels 9-24
    %tetrode_data(i,:) = amplifier_data(i,[trig:trig+intan_points-1]); %for data with electrodes on channel 1-16
end

tetrode_time = t_amplifier(1:intan_points); %time points for tetrode data

%create array of image time points
imaging_time = zeros(1,size(optical_data,2));
for i=1:size(optical_data,2)
    imaging_time(1,i) = i/framerate;
end


%Do bandwidth analysis

%Theta
lowF(1)=6;
highF(1)=14;

%Beta
lowF(2)=15;
highF(2)=30;

%Low gamma
lowF(3)=35;
highF(3)=55;

%High gamma
lowF(4)=65;
highF(4)=95;

%Names of bandwidths
bw_names{1}='Theta';
bw_names{2}='Beta';
bw_names{3}='Low gamma';
bw_names{4}='High gamma';

%%START DIEGOS CODE
% for elect_no=1:no_electrodes
%     decimated_LFP_logPtraces(elect_no,1)=mean(log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(1))&(LFPtime<=dFFtime(1)+(dt_dFF/2))),3);
%     decimated_LFP_logPtraces(elect_no,end)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(end)-(dt_dFF/2))&(LFPtime<=dFFtime(end))),3);
%     for ii_time=2:no_timepoints_dFF-1
%         decimated_LFP_logPtraces(:,ii_time)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(:,bwii,(LFPtime>dFFtime(ii_time)-(dt_dFF/2))&(LFPtime<=dFFtime(ii_time)+(dt_dFF/2))),3);
%     end
% end
%%END DIEGOS CODE

dec_n=100;
dectim=downsample(tetrode_time,dec_n);
%Setup the wavelet scales
%   scales = helperCWTTimeFreqVector(minfreq,maxfreq,f0,dt,NumVoices)
%   f0 - center frequency of the wavelet in cycles/unit time
%   dt - sampling interval
%   NumVoices - number of voices per octave

NumVoices=5;
minfreq=1;
dt_rhd=1/intan_freq;
f0=5/(2*pi);

a0 = 2^(1/NumVoices);
minscale = f0/(maxfreq*dt_rhd);
maxscale = f0/(minfreq*dt_rhd);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt_rhd;

%Now do the wavelet transform
for electNo=1:no_electrodes
    notch60HzFilt = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',intan_freq);
    this_LFP=zeros(1,size(tetrode_data(electNo,tetrode_time>0),2));
    this_LFP(1,:)=tetrode_data(electNo,tetrode_time>0);
    LFP=filtfilt(notch60HzFilt,this_LFP);


% dt_dFF=1/framerate;
% for elect_no=1:16
%     decimated_LFP_logPtraces(elect_no,1)=mean(log_P_timecourses_per_bw(elect_no,bwii,(tetrode_time>imaging_time(1))&(tetrode_time<=imaging_time(1)+(dt_dFF/2))),3);
%     decimated_LFP_logPtraces(elect_no,end)=mean(log_P_timecourses_per_bw(elect_no,bwii,(tetrode_time>imaging_time(end)-(dt_dFF/2))&(tetrode_time<=imaging_time(end))),3);
%     for ii_time=2:size(optical_data,2)-1
%         decimated_LFP_logPtraces(:,ii_time)=mean(log_P_timecourses_per_bw(:,bwii,(tetrode_time>imaging_time(ii_time)-(dt_dFF/2))&(tetrode_time<=imaging_time(ii_time)+(dt_dFF/2))),3);
%     end
% end
    
    decLFP=decimate(LFP,dec_n);
    decFs=intan_freq/dec_n;

    cwtLFP = cwtft({detrend(double(decLFP)),1/decFs},'wavelet','morl','scales',scales);
    Prev=abs(cwtLFP.cfs).^2;
    P=Prev(end:-1:1,:);

    frev=cwtLFP.frequencies;
    f=frev(end:-1:1);

    if electNo==1
        all_Power_timecourse=zeros(no_electrodes,length(f),length(dectim));
    end
    all_Power_timecourse(electNo,:,:)=P;

end

%Per electrode power plot
log_P_per_trial_timecourse_sub=zeros(length(f)*no_electrodes,length(dectim));
log_P_timecourses=zeros(no_electrodes,length(f),length(dectim));
log_P_timecourses_per_bw=zeros(no_electrodes,length(lowF),length(dectim));
mean_log_P_timecourse=zeros(length(f),length(dectim));
y_shift=0;
sy_shift=0;
shifted_freq=[];
for electNo=1:no_electrodes
    this_log_P_timecourse=zeros(length(f),length(dectim));
    this_log_P_timecourse(:,:)=10*log10(all_Power_timecourse(electNo,:,:));
    mean_log_P_timecourse(:,:)=mean_log_P_timecourse(:,:)+this_log_P_timecourse(:,:);
    log_P_per_trial_timecourse_sub(y_shift+1:y_shift+length(f),:)=this_log_P_timecourse;
    log_P_timecourses(electNo,:,:)=this_log_P_timecourse;
    for bwii=1:length(lowF)
        log_P_timecourses_per_bw(electNo,bwii,:)=mean(this_log_P_timecourse((f>=lowF(bwii))&(f<highF(bwii)),:),1);
    end
    shifted_freq(1,y_shift+1:y_shift+length(f))=f+(electNo-1)*f(end);
    y_shift=y_shift+length(f);
end
mean_log_P_timecourse=mean_log_P_timecourse/no_electrodes;

%     handles_per_file.file(handles_per_file.no_files).log_P_timecourses=log_P_timecourses;
log_P_timecourses_per_bw=log_P_timecourses_per_bw;
%     handles_per_file.file(handles_per_file.no_files).f=f;
log_P_time=dectim;

% %Plot the per-electrode timecourse
% figNo=0;
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% 
% maxLogPper=prctile(log_P_per_trial_timecourse_sub(:),99);
% minLogPper=prctile(log_P_per_trial_timecourse_sub(:),1);
% %Note: Diego added this on purpose to limit the range to 10 dB
% %This results in emphasizing changes in the top 10 dB
% if maxLogPper-minLogPper>12
%     minLogPper=maxLogPper-12;
% end
% 
% set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
% drg_pcolor(repmat(dectim,length(f)*no_electrodes,1)',repmat(shifted_freq,length(dectim),1),log_P_per_trial_timecourse_sub')
% 
% colormap jet
% shading interp
% %     caxis([minLogPper maxLogPper]);
% % xlim([0 576.8448])
% xlabel('Time (sec)', 'FontSize', 38)
% ylabel('Frequency', 'FontSize', 38);
% title('Tetrode Data: Power (dB, wavelet) Timecourse per Electrode', 'FontSize', 52)
% ax = gca;
% ax.FontSize = 28; 

% %Plot the mean wavelet power timecourse
% figNo=figNo+1;
% hFig=figure(figNo);
% 
% set(hFig, 'units','normalized','position',[.07 .1 .75 .3])
% drg_pcolor(repmat(dectim,length(f),1)',repmat(f,length(dectim),1),mean_log_P_timecourse')
% 
% colormap jet
% shading interp
% xlim([0 max(dectim)])
% xlabel('Time (s)', 'FontSize', 52)
% ylabel('Frequency (Hz)', 'FontSize', 52);
% title('Tetrode Data: Mean Power (dB, wavelet) Timecourse', 'FontSize', 64)
% ax = gca;
% ax.FontSize = 28; 
% c=colorbar('FontSize',28);
% c.Label.String = 'Power in dB';
 
seg_method=1; %If this variable is 0 the program uses random segmentation, otherwise it uses random circular shift

figNo=0;

k_fold=8;

% %Theta
% handles.lowF(1)=6;
% handles.highF(1)=14;
% 
% %Beta
% handles.lowF(2)=15;
% handles.highF(2)=30;
% 
% %Low gamma
% handles.lowF(3)=35;
% handles.highF(3)=55;
% 
% %High gamma
% handles.lowF(4)=65;
% handles.highF(4)=95;
% 
% %swr
% handles.lowF(4)=150;
% handles.highF(4)=250;
% 
% %Names of bandwidths
% handles.bw_names{1}='Theta';
% handles.bw_names{2}='Beta';
% handles.bw_names{3}='Low gamma';
% handles.bw_names{4}='High gamma';
% handles.bw_names{5}='swr';

%Ask user for the file with the traces
% [choiceFileName,choiceBatchPathName] = uigetfile({'*batch_per_file.mat'},'Select the .mat file for analysis');
% fprintf(1, ['\ndrgcmCrossPCArhddFF run for ' choiceFileName '\n\n']);

choiceBatchPathName = 'C:\Users\Connor\Desktop\'

% cd(choiceBatchPathName)
% load(choiceFileName)

%rho
ii_rho=0;
rho_dFF=[];
pval_rho_dFF=[];
rho_ddFF=[];
pval_rho_ddFF=[];
rho_fileNo=[];
rho_dFFtrace_no=[];
rho_electrode_no=[];
rho_bwii=[];

%Shuffled rho
ii_rhos=0;
rhos_dFF=[];
pval_rhos_dFF=[];
rhos_ddFF=[];
pval_rhos_ddFF=[];
rhos_fileNo=[];
rhos_dFFtrace_no=[];
rhos_electrode_no=[];
rhos_bwii=[];

%rhos for glm
ii_rho_glm=0;
rho_glm=[];
rho_glm_bwii=[];
rho_glm_elctrode_no=[];

%Now let's read Connor's data
optical_data = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\DFTemporal.csv'); %DF/F traces from CaImAn
highest_corr=zeros(4,16);
highest_corr_ROI=ones(4,16);
highest_rho_per_ROI=zeros(4,100);
elect_for_highest_rho_per_ROI=ones(4,100);
for fileNo=1:no_files
    for bwii=1:4
        no_timepoints_dFF=size(optical_data,2);
        %         if use_raw==1
        %             dFFtraces=handles_per_file.file(fileNo).dFFtraces;
        %         else
        %             dFFtraces=handles_per_file.file(fileNo).dFFtraces_inferred;
        %         end
        dFFtraces=optical_data;
        no_traces_dFF=size(optical_data,1);
        
        LFPtime=tetrode_time;
        dFFtime=imaging_time;
        dFFtime=dFFtime(1:size(dFFtraces,2));
        no_timepoints_dFF=size(dFFtraces,2);
        %decimated_LFP_logPtraces=zeros(handles_per_file.file(fileNo).no_electrodes,no_timepoints_dFF);
        dt_dFF=1/framerate;
%         for elect_no=1:handles_per_file.file(fileNo).no_electrodes
%             decimated_LFP_logPtraces(elect_no,1)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(1))&(LFPtime<=dFFtime(1)+(dt_dFF/2))),3);
%             decimated_LFP_logPtraces(elect_no,end)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(elect_no,bwii,(LFPtime>dFFtime(end)-(dt_dFF/2))&(LFPtime<=dFFtime(end))),3);
%             for ii_time=2:no_timepoints_dFF-1
%                 decimated_LFP_logPtraces(:,ii_time)=mean(handles_per_file.file(fileNo).log_P_timecourses_per_bw(:,bwii,(LFPtime>dFFtime(ii_time)-(dt_dFF/2))&(LFPtime<=dFFtime(ii_time)+(dt_dFF/2))),3);
%             end
%         end

        decimated_LFP_logPtraces = decLFP;
        optical_interp = zeros(size(optical_data,1), size(log_P_timecourses_per_bw,3));
        fprintf('Beginning to interpolate optical data \n')
        for i=1:size(optical_data,1)
            optical_interp(i,:)=interp1(imaging_time, optical_data(i,:), dectim, 'linear','extrap');
        end
        fprintf('Done interpolating optical data, moving on to correlation analysis \n')
        
        %Save the decimated data
%         save([choiceFileName(1:end-19) '_dFF.txt'],'dFFtraces','-ascii')
%         save([choiceFileName(1:end-19) '_LFP_' handles.bw_names{bwii} '.txt'],'decimated_LFP_logPtraces','-ascii')
%         dFFtraces_t=dFFtraces';
%         decimated_LFP_logPtraces_t=decimated_LFP_logPtraces';
%         save([choiceFileName(1:end-19) '_dec_' handles.bw_names{bwii} '.mat'],'dFFtraces_t','decimated_LFP_logPtraces_t','dFFtime')
        
        %Calculate rho
        mask_dFF=~isnan(decimated_LFP_logPtraces(1,:));
        for elect_no=1:no_electrodes
            for trace_no=1:no_traces_dFF
                ii_rho=ii_rho+1;
                [rho_temp,pval_rho_dFF_temp] = corrcoef(optical_interp(trace_no,:),log_P_timecourses_per_bw(elect_no,bwii,:));
                rho_dFF(ii_rho) = rho_temp(1,2);
                pval_rho_dFF(ii_rho) = pval_rho_dFF_temp(1,2);
                rho_fileNo(ii_rho)=fileNo;
                rho_dFFtrace_no(ii_rho)=trace_no;
                rho_electrode_no(ii_rho)=elect_no;
                rho_bwii(ii_rho)=bwii;
                
                if highest_corr(bwii,elect_no)<rho_dFF(ii_rho)
                    highest_corr(bwii,elect_no)=rho_dFF(ii_rho);
                    highest_corr_ROI(bwii,elect_no)=trace_no;
                end
                
                if highest_rho_per_ROI(bwii,trace_no)<rho_dFF(ii_rho)
                    highest_rho_per_ROI(bwii,trace_no)=rho_dFF(ii_rho);
                    elect_for_highest_rho_per_ROI(bwii,trace_no)=elect_no;
                end
                
            end
        end
        
        %Shuffle dFF
        no_shuffles=100;
        if seg_method==0
            %Now calculate rho for shuffled dFF by random segmentation
            no_segments=10;
            
            segments=[1:no_segments];
            per_seg=perms(segments);
            seg_ii=randi(size(per_seg,1),1,10000);
            ii_used=0;
            seg_length=floor(sum(mask_dFF)/no_segments);
            
            for shuf_ii=1:no_shuffles
                found_shuffle=0;
                while found_shuffle==0
                    ii_used=ii_used+1;
                    if sum(per_seg(seg_ii(ii_used),:)==segments)==0
                        found_shuffle=1;
                    end
                end
                
                for elect_no=1:no_electrodes
                    for trace_no=1:no_traces_dFF
                        shdFFtrace=zeros(1,no_segments*seg_length);
                        this_per_ii=per_seg(seg_ii(ii_used),:);
                        for ii=1:no_segments
                            shdFFtrace(1,(ii-1)*seg_length+1:ii*seg_length)=dFFtraces(trace_no,(this_per_ii(ii)-1)*seg_length+1:this_per_ii(ii)*seg_length);
                        end
                        ii_rhos=ii_rhos+1;
                        [rhos_dFF_temp,pval_rhos_dFF_temp] = corrcoef(shdFFtrace(1,:)',log_P_timecourses_per_bw(elect_no,bwii,1:length(shdFFtrace))');
                        rhos_dFF(ii_rhos)=rhos_dFF_temp(1,2);
                        pval_rhos_dFF(ii_rhos)=pval_rhos_dFF_temp(1,2);
                        rhos_fileNo(ii_rhos)=fileNo;
                        rhos_dFFtrace_no(ii_rhos)=trace_no;
                        rhos_electrode_no(ii_rhos)=elect_no;
                        rhos_bwii(ii_rhos)=bwii;
                    end
                end
            end
            
        else
            %Now calculate rho for shuffled dFF using a circular shift
            
            
            shift_ii=ceil(0.1*size(dFFtraces,2)) + randi(ceil(0.9*size(dFFtraces,2)),1,1000);
            no_timepoints=size(dFFtraces,2);
            
            for shuf_ii=1:no_shuffles
                for elect_no=1:no_electrodes
                    for trace_no=1:no_traces_dFF
                        shdFFtrace=zeros(1,no_timepoints);
                        
                        %Shift forward
                        shdFFtrace(1,shift_ii(shuf_ii):end)=dFFtraces(trace_no,1:no_timepoints-shift_ii(shuf_ii)+1);
                        shdFFtrace(1,1:shift_ii(shuf_ii)-1)=dFFtraces(trace_no,no_timepoints-shift_ii(shuf_ii)+2:end);
                        
                        ii_rhos=ii_rhos+1;
                        [rhos_dFF_temp,pval_rhos_dFF_temp] = corrcoef(shdFFtrace(1,:)',log_P_timecourses_per_bw(elect_no,bwii,1:length(shdFFtrace)));
                        rhos_dFF(ii_rhos)=rhos_dFF_temp(1,2);
                        pval_rhos_dFF(ii_rhos)=pval_rhos_dFF_temp(1,2);
                        rhos_fileNo(ii_rhos)=fileNo;
                        rhos_dFFtrace_no(ii_rhos)=trace_no;
                        rhos_electrode_no(ii_rhos)=elect_no;
                        rhos_bwii(ii_rhos)=bwii;
                    end
                end
            end
        end
        %          %Fit LFP with dFFtraces
        %         order=3;
        %         framelen=31;
        %         mask_dFF=~isnan(decimated_LFP_logPtraces(1,:));
        %         for elect_no=1:handles_per_file.file(fileNo).no_electrodes
        %             this_decimated_LFP_logPtraces=zeros(sum(mask_dFF),1);
        %             this_decimated_LFP_logPtraces(:,1)=decimated_LFP_logPtraces(elect_no,mask_dFF)';
        %             sgf_this_decimated_LFP_logPtraces = sgolayfilt(this_decimated_LFP_logPtraces,order,framelen);
        %             %Now do a k-fold glm fit where the data left out (1/8th) is predicted from
        %             %the data that was used for the glm (the other 7/8ths)
        %             k_fold_chunk=floor(length(sgf_this_decimated_LFP_logPtraces)/k_fold);
        %             LFPlog_pred=zeros(length(sgf_this_decimated_LFP_logPtraces),1);
        %             LFPlog_CI=zeros(length(sgf_this_decimated_LFP_logPtraces),2);
        %             for k_fold_ii=1:k_fold
        %                 mask_test_dFFtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 mask_test_dFFtraces((k_fold_ii-1)*k_fold_chunk+1:k_fold_ii*k_fold_chunk)=1;
        %                 mask_test_dFFtraces=logical(mask_test_dFFtraces);
        %                 mdl=fitglm(dFFtraces(:,~mask_test_dFFtraces)',sgf_this_decimated_LFP_logPtraces(~mask_test_dFFtraces),'linear');
        %                 [LFPlog_pred(mask_test_dFFtraces),LFPlog_CI(mask_test_dFFtraces,:)] = predict(mdl,dFFtraces(:,mask_test_dFFtraces)');
        %             end
        %             ii_rho_glm=ii_rho_glm+1;
        %             rho_glm_elctrode_no(ii_rho_glm)=elect_no;
        %             rho_glm_bwii(ii_rho_glm)=bwii;
        %             rho_glm(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
        %
        %
        %             if figNo<5
        %                 figNo=figNo+1;
        %                 try
        %                     close(figNo)
        %                 catch
        %                 end
        %                 hFig=figure(figNo);
        %                 hold on
        %                 set(hFig, 'units','normalized','position',[.3 .3 .6 .3])
        %                 plot(dFFtime,sgf_this_decimated_LFP_logPtraces,'-k','LineWidth',3)
        %                 plot(dFFtime+3,LFPlog_pred,'-r')
        %                 title(['glm fit of logLFP of electrode No ' num2str(elect_no) ' for ' handles.bw_names{bwii}])
        %                 fprintf(1, ['\nProcessed glm for ' handles.bw_names{bwii} ' electrode No ' num2str(elect_no) '\n']);
        %             end
        %
        %             %Do shuffled glm
        %             LFPlog_pred=zeros(length(sgf_this_decimated_LFP_logPtraces),1);
        %             LFPlog_CI=zeros(length(sgf_this_decimated_LFP_logPtraces),2);
        %
        %             for k_fold_ii=1:k_fold
        %                 mask_test_dFFtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 mask_test_dFFtraces((k_fold_ii-1)*k_fold_chunk+1:k_fold_ii*k_fold_chunk)=1;
        %                 mask_test_dFFtraces=logical(mask_test_dFFtraces);
        %                 perm_k_fold_ii=randperm(k_fold);
        %                 sh_sgf_this_decimated_LFP_logPtraces=zeros(1,length(sgf_this_decimated_LFP_logPtraces));
        %                 for jj=1:k_fold
        %                     sh_sgf_this_decimated_LFP_logPtraces((jj-1)*k_fold_chunk+1:jj*k_fold_chunk)=...
        %                         sgf_this_decimated_LFP_logPtraces((perm_k_fold_ii(jj)-1)*k_fold_chunk+1:perm_k_fold_ii(jj)*k_fold_chunk);
        %                 end
        %                 mdl=fitglm(dFFtraces(:,~mask_test_dFFtraces)',sh_sgf_this_decimated_LFP_logPtraces(~mask_test_dFFtraces),'linear');
        %                 [LFPlog_pred(mask_test_dFFtraces),LFPlog_CI(mask_test_dFFtraces,:)] = predict(mdl,dFFtraces(:,mask_test_dFFtraces)');
        %             end
        %
        %             rho_glm_sh(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
        %
        %         end
        
        %         %Calculate PCs per timepoint
        %         for ii_time=1:no_timepoints_dFF
        %             these_logPs=zeros(1,handles_per_file.file(fileNo).no_electrodes);
        %             these_logPs(1,:)=decimated_LFP_logPtraces(:,ii_time);
        %             this_PCs_logP=[];
        %             [coeff_logP,this_PCs_logP,latent_logP]=pca(these_logPs');
        %             if ii_time==1
        %                 PCs_logP=zeros(length(this_PCs_logP),no_timepoints_dFF);
        %             end
        %             PCs_logP(:,ii_time)=this_PCs_logP;
        %
        %             these_dFFs=zeros(1,no_traces_dFF);
        %             these_dFFs(1,:)=dFFtraces(:,ii_time);
        %             this_PCs_dFF=[];
        %             [coeff_logP,this_PCs_dFF,latent_logP]=pca(these_dFFs');
        %              if ii_time==1
        %                 PCs_dFF=zeros(length(this_PCs_dFF),no_timepoints_dFF);
        %             end
        %             PCs_dFF(:,ii_time)=this_PCs_dFF;
        %         end
        %
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %
        %         PC1_dFF=PCs_dFF(1,:)';
        %         PC1_logP=PCs_logP(1,:)';
        %         PC_mask=(~isnan(PC1_dFF))&(~isnan(PC1_logP));
        %         PC1_dFF=PC1_dFF(PC_mask);
        %         PC1_logP=PC1_logP(PC_mask);
        %
        %         plot(PC1_dFF,PC1_logP,'.k')
        %         xlabel('PC1 for dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(PC1_dFF,PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'PC1 dFF vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
        %Now do dFF derivative
        %Convolve lick_freq using a window of 0.9 sec
        %         no_conv_points=5;
        %         conv_win=ones(1,no_conv_points);
        %         conv_dFFtraces=zeros(handles_per_file.file(fileNo).no_electrodes,no_timepoints_dFF);
        %
        %         for trace_no=1:no_traces_dFF
        %             conv_dFFtraces(trace_no,:)=conv(dFFtraces(trace_no,:),conv_win,'same')/no_conv_points;
        %         end
        %
        %         ddFFtraces=zeros(no_traces_dFF,size(dFFtraces,2));
        %         ddFFtraces(:,2:end)=(conv_dFFtraces(:,2:end)-conv_dFFtraces(:,1:end-1))/dt_dFF;
        %         ddFFtraces(:,1)=ddFFtraces(:,2);
        %         pffft=1;
        %               %Calculate PCs per timepoint
        %         for ii_time=1:no_timepoints_dFF
        %
        %
        %             these_ddFFtraces=zeros(1,no_traces_dFF);
        %             these_ddFFtraces(1,:)=ddFFtraces(:,ii_time);
        %             this_PCs_ddFF=[];
        %             [coeff_logP,this_PCs_ddFF,latent_logP]=pca(these_ddFFtraces');
        %              if ii_time==1
        %                 PCs_ddFF=zeros(length(this_PCs_ddFF),no_timepoints_dFF);
        %             end
        %             PCs_ddFF(:,ii_time)=this_PCs_ddFF;
        %         end
        %
        %         PCs_ddFF=PCs_ddFF(:,PC_mask);
        %         PC1_ddFF=PCs_ddFF(1,:);
        %
        %         figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %         plot(PC1_ddFF,PC1_logP,'.k')
        %         xlabel('PC1 for derivative of dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(PC1_dFF,PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'PC1 dFF derivative vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
        %Now add all derivatives
        
        %         sumddFFtraces=sum(ddFFtraces,1);
        %          sumddFFtraces=sumddFFtraces(1,PC_mask);
        %          figNo=figNo+1;
        %         try
        %             close(figNo)
        %         catch
        %         end
        %         hFig=figure(figNo);
        %         set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
        %
        %         hold on
        %
        %
        %         plot(sumddFFtraces,PC1_logP,'.k')
        %         xlabel('Sum of derivatives fot dFF')
        %         ylabel('PC1 for logP')
        %         title(['PC1 plot for ' handles.bw_names{bwii}])
        %
        %         [rho,pval] = corr(sumddFFtraces',PC1_logP);
        %
        %         fprintf(1, ['\nFor file No %d ' handles.bw_names{bwii} 'sum of dFF derivatives vs PC1 logP rho= %d, p value= %d\n'],fileNo,rho, pval);
        %
    end
end
%
% %Plot correlation between glm prediction and LFP vs shuffled
% edges=[-1:0.05:1];
% rand_offset=0.8;
%
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% hFig=figure(figNo);
% set(hFig, 'units','normalized','position',[.3 .3 .6 .3])
% hold on
% bar_offset=0;
%
% glm_rho_glm=[];
% glm_ii=0;
%
% id_ii=0;
% input_data=[];
%
% for bwii=1:4
%
%     rho_glm_elctrode_no(ii_rho_glm)=elect_no;
%     rho_glm_bwii(ii_rho_glm)=bwii;
%     rho_glm(ii_rho_glm)=corr(LFPlog_pred,sgf_this_decimated_LFP_logPtraces);
%
%
%     %Shuffled rho
%     bar_offset=bar_offset+1
%     bar(bar_offset,mean(rho_glm_sh(rho_glm_bwii==bwii)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0.7 0.7 0.7])
%
%     %Violin plot
%     these_rhos=rho_glm_sh(rho_glm_bwii==bwii);
%     [mean_out, CIout]=drgViolinPoint(these_rhos...
%         ,edges,bar_offset,rand_offset,'k','k',3);
%
%     %Shuffled
%     glm_rho_glm.data(glm_ii+1:glm_ii+length(these_rhos))=these_rhos;
%     glm_rho_glm.bwii(glm_ii+1:glm_ii+length(these_rhos))=bwii*ones(1,length(these_rhos));
%     glm_rho_glm.shuffled(glm_ii+1:glm_ii+length(these_rhos))=ones(1,length(these_rhos));
%
%     glm_ii=glm_ii+length(these_rhos);
%
%     id_ii=id_ii+1;
%     input_data(id_ii).data=these_rhos;
%     input_data(id_ii).description=['rho for shuffled ' handles.bw_names{bwii}];
%
%
%     %rho
%     bar_offset=bar_offset+1
%     bar(bar_offset,mean(rho_glm(rho_glm_bwii==bwii)),'LineWidth', 3,'EdgeColor','none','FaceColor',[0 0.7 0.7])
%
%     %Violin plot
%     these_rhos=rho_glm(rho_glm_bwii==bwii);
%     [mean_out, CIout]=drgViolinPoint(these_rhos...
%         ,edges,bar_offset,rand_offset,'k','k',3);
%
%
%     bar_offset=bar_offset+1
%
%     %Original data
%     glm_rho_glm.data(glm_ii+1:glm_ii+length(these_rhos))=these_rhos;
%     glm_rho_glm.bwii(glm_ii+1:glm_ii+length(these_rhos))=bwii*ones(1,length(these_rhos));
%     glm_rho_glm.shuffled(glm_ii+1:glm_ii+length(these_rhos))=zeros(1,length(these_rhos));
%
%     glm_ii=glm_ii+length(these_rhos);
%
%     id_ii=id_ii+1;
%     input_data(id_ii).data=these_rhos;
%     input_data(id_ii).description=['rho for original ' handles.bw_names{bwii}];
%
%
% end
%
% title('Correlation between glm predict and LFP power')
% ylabel('rho')
%
% xticks([1.5 4.5 7.5 10.5])
% xticklabels({handles.bw_names{1}, handles.bw_names{2}, handles.bw_names{3}, handles.bw_names{4}})
%
% %Perform the glm
% fprintf(1, ['glm for rho between LFP log power and predicted using dFF\n'])
%
%
% tbl = table(glm_rho_glm.data',glm_rho_glm.bwii',glm_rho_glm.shuffled',...
%     'VariableNames',{'rho','bwii','shuffled'});
% mdl = fitglm(tbl,'rho~bwii+shuffled+bwii*shuffled'...
%     ,'CategoricalVars',[2,3])
%
%
% %Do the ranksum/t-test
% fprintf(1, ['\n\nRanksum or t-test p values for rho between LFP log power and predicted using dFF\n'])
%
% try
%     [output_data] = drgMutiRanksumorTtest(input_data);
% catch
% end

%Plot correlation histograms
edges=[-0.4:0.01:0.4];

pFDR=drsFDRpval(pval_rho_dFF);
fprintf(1,'\n\npFDR for rho p value  = %d\n',pFDR);

pFDRs=drsFDRpval(pval_rhos_dFF);
fprintf(1, ['\n\npFDR for shuffled rho p value = %d\n'],pFDR);

for bwii=1:4
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .6])
    edges=[-0.31:0.01:0.31];
    %subplot(2,1,1)
    hold on
    histogram(rho_dFF(rho_bwii==bwii),edges);
    %histogram(rho_dFF(rho_bwii==bwii));

    title('Original')
    xlabel('rho')
    
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .6])
    
    %Plot shuffled correlation histograms
    
   
    

    hold on
    edges=[-0.31:0.01:0.31];
    histogram(rhos_dFF(rhos_bwii==bwii),edges,'FaceColor',[0.3010 0.7450 0.9330]);
    %histogram(rhos_dFF(rhos_bwii==bwii),'FaceColor',[0.3010 0.7450 0.9330]);

    title('Shuffled')
    xlabel('rho')
    
    %sgtitle(['rho for dFF x LFP log P for ' bw_names{bwii} ])
    
    [h,p] = kstest2(rhos_dFF(rhos_bwii==bwii),rho_dFF(rho_bwii==bwii),'Alpha',0.01);
    fprintf(1, ['\np value for thr KS test difference between original and shuffled ' bw_names{bwii}  ' = %d\n\n'],p);
    
    
    figNo=figNo+1;
    try
        close(figNo)
    catch
    end
    hFig=figure(figNo);
    set(hFig, 'units','normalized','position',[.3 .3 .3 .3])
    hold on
    
    [f_rho,x_rho] = drg_ecdf(rho_dFF);
    plot(x_rho,f_rho,'b','LineWidth',3)
    
    [f_rhos,x_rhos] = drg_ecdf(rhos_dFF);
    plot(x_rhos,f_rhos,'Color',[0.3010 0.7450 0.9330],'LineWidth',3)
    
    text(0.2,0.6,'Original','Color','b','FontSize',12)
    text(0.2,0.5,'Shuffled','Color',[0.3010 0.7450 0.9330],'FontSize',12)

    title(['Cumulative probability for rho ' bw_names{bwii}])
    xlabel('rho')
    ylabel('Cumulative probability')
    
end


figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.3 .3 .4 .4])

hold on

ii=0;
for bwii=1:4
    ii=ii+1;
    bar(ii,100*sum(pval_rho_dFF(rho_bwii==bwii)<pFDR)/sum(rho_bwii==bwii),'r')
    ii=ii+1;
    bar(ii,100*sum(pval_rhos_dFF(rhos_bwii==bwii)<pFDRs)/sum(rhos_bwii==bwii),'b')
    ii=ii+1;
end

xticks([1.5 4.5 7.5 10.5])
xticklabels({'Theta','Beta','Low gamma','High gamma'})
ylim([0 110])
ylabel('Percent significant rho')
title('Percent significant correlations')
text(10,38,'Original','Color','red','FontSize',12)
text(10,36,'Shuffled','Color','blue','FontSize',12)


figNo=figNo+1;
try
    close(figNo)
catch
end
hFig=figure(figNo);
set(hFig, 'units','normalized','position',[.3 .3 .3 .3])

hold on

ii=0;
for bwii=1:4
    bar(ii,100*sum(pval_rho_dFF(rho_bwii==bwii)<pFDR)/sum(rho_bwii==bwii)-100*sum(pval_rhos_dFF(rhos_bwii==bwii)<pFDRs)/sum(rhos_bwii==bwii),'b')
    ii=ii+2;
end

xticks([0 1 2 3])
xticklabels({'Theta','Beta','Low gamma','High gamma'})
ylim([0 110])
ylabel('Percent significant rho')
title('Percent significant correlations')


writematrix(highest_corr,[choiceBatchPathName 'highest_corr.csv'])
writematrix(highest_corr_ROI,[choiceBatchPathName 'highest_corr_ROI.csv'])
writematrix(highest_rho_per_ROI',[choiceBatchPathName 'highest_rho_per_ROI.csv'])
writematrix(elect_for_highest_rho_per_ROI',[choiceBatchPathName 'elect_for_highest_rho_per_ROI.csv'])


pffft=1;
