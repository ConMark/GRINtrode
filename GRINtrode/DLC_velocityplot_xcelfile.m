function DLC_velocityplot_xcelfile

data = readmatrix('E:\TL12_rencodeDLC_resnet50_Onecam_WholeMouseMar24shuffle1_500000.csv');

xpos = data(:,8);
ypos = data(:,9);
prob = data(:,10);
framerate = 30; %framerate of behavioral video (frames/second)

velocity = zeros(size(xpos,1),1);
velocity(1,1)=0;

time = zeros(size(xpos,1),1);
time(1,1)=0;

for i = 2:size(xpos,1)
    if prob(i,1) > 0.9
        velocity(i,1) = sqrt((xpos(i)-xpos(i-1))^2+(ypos(i)-ypos(i-1))^2)*framerate;
    else
        velocity(i,1) = NaN;
    end
    time(i,1) = i/framerate;
end

figNo=0;
figNo=figNo+1;
figure(figNo)
plot(time, velocity)
size(velocity)

trimmed_velocity = velocity(75:14586,1);
num_frames=size(trimmed_velocity,1);
for i = 1:num_frames
    trimmed_time(i,1) = i/framerate;
end



n = 15; % average every n values

for i = 1:size(trimmed_velocity,1)/n
    averaged_time(i,1) = (n/framerate)*i;
end

averaged_velocity = arrayfun(@(i) mean(trimmed_velocity(i:i+n-1)),1:n:length(trimmed_velocity)-n+1)'; % the averaged vector
figNo=figNo+1;
figure(figNo)
plot(averaged_time,averaged_velocity)
size(averaged_velocity)
ax = gca;
ax.FontSize = 54;
xlabel('Time (seconds)', 'FontSize', 60);
ylabel('Velocity (pixels/second)', 'FontSize', 60);
xlim([0 max(averaged_time)])
max(averaged_time)