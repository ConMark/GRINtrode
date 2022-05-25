function DLC_velocityplot_matfile

data = load('E:\T12_pos_vel_acc_mat.mat');

xpos = data.data.mhbcenter(:,1);
ypos = data.data.mhbcenter(:,2);
framerate = 30; %frames/second

velocity = zeros(size(xpos,1),1);
velocity(1,1)=0;

time = zeros(size(xpos,1),1);
time(1,1)=0;

for i = 2:size(xpos,1)
    velocity(i,1) = sqrt((xpos(i)-xpos(i-1))^2+(ypos(i)-ypos(i-1))^2)*framerate;
    time(i,1) = i/framerate;
end

figNo=0;
figNo=figNo+1;
figure(figNo)
plot(time, velocity)


