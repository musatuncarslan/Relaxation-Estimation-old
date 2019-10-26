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

fs = MPIparams.fs; f_d = MPIparams.f_drive; 
t = (0:fs*5/f_d-1)/fs; % time axis for error calculation


slewRate = 0:0.5:20; % T/s
Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
G=MPIparams.Gzz/Physicsparams.mu0; % gradient
MPIparams.driveMag=Hp/G; % extent of the drive field

R = slewRate/MPIparams.Gzz;


arg = 2*pi*f_d;
p = zeros(1, 3);
for k=0:length(p)-1
    p(k+1) = (-1)^k/factorial(2*k+1)*arg^(2*k+1);
end
p = flip(p)*MPIparams.driveMag;
p = upsample(p, 2);

t_idx = 51;
for R_idx = 1:length(R)
    p_temp = p;
    p_temp(end-1) = p_temp(end-1) - R(R_idx);
    roots_p = roots(p_temp);
    real_idx = imag(roots_p) == 0;
    real_roots = roots_p(real_idx);

    del_t(R_idx) = real_roots(end);
    err(R_idx) = abs((-MPIparams.driveMag*cos(2*pi*f_d*t(t_idx)) + R(R_idx)*t(t_idx)) - ...
        (-MPIparams.driveMag*cos(2*pi*f_d*(t(t_idx)+del_t(R_idx))) + R(R_idx)*(t(t_idx) + del_t(R_idx))));  
end

figure; 
subplot(2,1,1); plot(slewRate, del_t*1000); axis tight;
xlabel('Slew Rate (T/s)'); ylabel('\Delta t (ms)');

subplot(2,1,2); plot(slewRate, err); axis tight;
xlabel('Slew Rate (T/s)'); ylabel('Absolute Error');    
% ylim([0, 10^-3])

figure; 
plot(t, -MPIparams.driveMag*cos(2*pi*f_d*t)+R(end)*t)
hold on
plot(t(t_idx), -MPIparams.driveMag*cos(2*pi*f_d*t(t_idx)) + R(end)*t(t_idx), '*')
plot(t(t_idx)+del_t(end), -MPIparams.driveMag*cos(2*pi*f_d*(t(t_idx)+del_t(end))) + R(end)*(t(t_idx) + del_t(end)), '*')
xlabel('Time (s)'); ylabel('FFP Position (m)');  