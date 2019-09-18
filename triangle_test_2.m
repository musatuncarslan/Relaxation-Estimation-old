clear all
close all
clc

format long

FOV_x = 1;
FOV_z = 0.6;
time = 1;
numPeriod = 4*2;
fs = 4e4; p = 1/numPeriod; 


traversedFOVz = [-0.03 0.03]; % traversed fov in z-axis (m)
traversedFOVx = [-0.01 0.01]; % traversed fov in x-axis (m)

robotSpeed_z = FOV_z/time; % robot arm movement speed in z direction (m/s)
bFOVz = -FOV_z/2; % beginning point of FOV in z-axis (m)
t1_z = round((traversedFOVz(1)-bFOVz)/robotSpeed_z, 5);
t2_z = round((traversedFOVz(2)-bFOVz)/robotSpeed_z, 5);

t1_x = [-p*(0.5-traversedFOVx(1)/FOV_x) p*(0.5-traversedFOVx(1)/FOV_x)]+p;
t2_x = [-p*(0.5-traversedFOVx(2)/FOV_x) p*(0.5-traversedFOVx(2)/FOV_x)]+p;

t = -p:1/(fs):p-1/(fs);
x = (1-abs(t/p))*FOV_x - FOV_x/2;

x_c = [];
t_c = [];
for k=1:numPeriod*time/2
    x_c = [x_c x];
    t_c = [t_c t+(k-1)*p*2];
    t1_x(k, :) = [-p*(0.5-traversedFOVx(1)/FOV_x) p*(0.5-traversedFOVx(1)/FOV_x)]+(k-1)*p*2;
    t2_x(k, :) = [-p*(0.5-traversedFOVx(2)/FOV_x) p*(0.5-traversedFOVx(2)/FOV_x)]+(k-1)*p*2;
end
t_c = t_c+p;
t1_x = round(t1_x + p, 5)
t2_x = round(t2_x + p, 5)

for k=1:size(t1_x, 1)
    
end


t_2 = (0:time*fs-1)/fs;
x_2 = FOV_x*(2*abs(2*(t_2/(p*2) - floor(t_2/(p*2) + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time
z = FOV_z/time*t_2 - FOV_z/2; % robot arm movement in z direction w.r.t. time

figure; plot(t_c, x_c, 'linewidth', 2);
% figure; plot(t_c, x_c, 'linewidth', 2); hold on; plot(t_2, x_2, 'linewidth', 2); legend('triangular wave', 'triangular function')
% figure; plot(x_c-x_2); legend('difference')
% 
% figure; plot(x_c, z, '*', 'linewidth', 2); hold on; plot(x_2, z);
