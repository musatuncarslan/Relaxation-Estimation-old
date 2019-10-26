clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 20); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams, 256, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

slewRate = 0:0.25:20; % T/s
Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
G=MPIparams.Gzz/Physicsparams.mu0; % gradient
MPIparams.driveMag=Hp/G; % extent of the drive field

R = slewRate/MPIparams.Gzz;


f_d = MPIparams.f_drive; 
fs = 100e6;
t = (0:fs*1/f_d-1)/fs; % time axis for error calculation



for l=1:length(R)
    del_t(l) = findShift(R(l), f_d, MPIparams);
    err(l) = abs(MPIparams.driveMag*cos(2*pi*f_d/4/f_d) - (MPIparams.driveMag*cos(2*pi*f_d*(1/(4*f_d) + del_t(l))) + R(l)*del_t(l)));

end
figure; 
subplot(2,1,1); plot(slewRate, del_t); title('\Delta t w.r.t. Slew Rate'); xlabel('R_s (T/s)'); ylabel('\Delta t (s)')
subplot(2,1,2); plot(slewRate, err); title('Error of Analytic solution'); xlabel('R_s (T/s)'); ylabel('Absolute Error')

figure;
idx = length(R);
plot(t, MPIparams.driveMag*cos(2*pi*f_d*t) + R(idx)*t); hold on;
plot(1/4/f_d, MPIparams.driveMag*cos(2*pi*f_d/4/f_d) + R(idx)/4/f_d, '*');
plot(1/4/f_d+del_t(idx), MPIparams.driveMag*cos(2*pi*f_d*(1/4/f_d+del_t(idx))) + R(idx)*(1/4/f_d+del_t(idx)), '*');
title(['FFP Position w.r.t. time for R_s = ' num2str(slewRate(idx)) ' T/s'])
ylabel('z-axis position (m)'); xlabel('t (s)'); legend('FFP movement', 'Position @ max speed', 'Analytic sol. to pos. @ max speed')
