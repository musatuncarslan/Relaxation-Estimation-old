clear all
close all
clc

fs = 0.5e6;
FOV_x = 0.05;
FOV_z = 0.06;
time = 3; % time to traverse whole FOV

traversedFOVz = [-0.03 0.03]; % traversed fov in z-axis (m)
traversedFOVx = [-0.005 0]; % traversed fov in x-axis (m)

robotSpeed_z = FOV_z/time; % robot arm movement speed in z direction (m/s)
bFOVz = -FOV_z/2; % beginning point of FOV in z-axis (m)
t1_z = round((traversedFOVz(1)-bFOVz)/robotSpeed_z, 2);
t2_z = round((traversedFOVz(2)-bFOVz)/robotSpeed_z, 2);

numPeriod = 5;
% robotSpeed_x = (FOV_x*numPeriod*time)/time; % robot arm movement speed in z direction (m/s)
% bFOVx = -FOV_x/2; % beginning point of FOV in z-axis (m)
% t1_z = round((traversedFOVx(1)-bFOVx)/robotSpeed_x, 2);
% t2_z = round((traversedFOVx(2)-bFOVx)/robotSpeed_x, 2);


total_time = round(t2_z-t1_z, 2);



p = 1/numPeriod;

t = (0:total_time*fs-1)/fs;

x = FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time
z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time

plot(t, x); hold on; plot(t, z);

figure; plot(x, z); xlabel('x-axis'); ylabel('z-axis')
xlim([-0.025 0.025]); ylim([-0.03 0.03])