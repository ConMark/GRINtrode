function align_spikes_and_rhd
%time syncs Spikes.csv from CaImAn and raw tetrode data from Intan

%run read_Intan_RHD2000_file.m first

%imports RHD data
T = evalin('base','whos');
for ii = 1:length(T)
   C_ =  evalin('base',[T(ii).name ';']);
   eval([T(ii).name,'=C_;']);
end
clear T C_ ii

intan_freq = frequency_parameters.amplifier_sample_rate %intan sampling rate in Hz
no_electrodes=16;

optical_spikes = readmatrix('C:\Users\Connor\Documents\Work\Paper\157565_9_27_TL2\Spikes.csv'); %DF/F traces from CaImAn
framerate = 1/0.25984; %framerate in frames/second

trig=find(board_dig_in_data(6,:),1,'first'); %find trigger input from imaging start
imaging_start_time = trig*(1/intan_freq) %when imaging started after starting intan recording
imaging_length_seconds = size(optical_spikes,2)/framerate %length of timelapse
intan_points = int64(imaging_length_seconds*20000) %number of intan data points

%create tetrode data array
tetrode_data = zeros(no_electrodes,intan_points);
for i = 1:no_electrodes
    tetrode_data(i,:) = amplifier_data(i+8,[trig:trig+intan_points-1]); %for data with electrodes on channels 9-24
    %tetrode_data(i,:) = amplifier_data(i,[trig:trig+intan_points-1]); %for data with electrodes on channel 1-16
end

tetrode_time = t_amplifier(1:intan_points); %time points for tetrode data

%create array of image time points
imaging_time = zeros(1,size(optical_spikes,2));
for i=1:size(optical_spikes,2)
    imaging_time(1,i) = i/framerate;
end

figure
subplot(1,2,1)
plot(tetrode_time,tetrode_data(1,:))
subplot(1,2,2)
plot(imaging_time,optical_spikes(1,:))