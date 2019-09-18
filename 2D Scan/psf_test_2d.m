clear all
close all
clc

format long;

% gradients (T/m)
Gxx = -4;
Gyy = 0;
Gzz = 4;
diameter = 30;

Bp = 15e-3; % Drive field (T)
f_drive = 10e3; % drive field frequency
mu0=1.256637*10^-6; % permaebility of vacuum
Hp=Bp/mu0; % magnetization moment
G=Gzz/mu0; % gradient
driveMag=Hp/G; % extent of the drive field

fs = 2e6; % sample frequency of the physical world (Hz)
fs_mpi = 2e6; % sample frequency of the MPI system (Hz)
FOV_z = 0.06; % FOV in z-axis (meters) (bore axis)
FOV_x = 0.05; % FOV in x-axis (meters)
time = 2; % time (seconds)

% generate FFP movements
[FFP_x, FFP_z, FFP_speed, FFP_angle, x, z] = FFPtrajectory(FOV_x, FOV_z, fs, f_drive, time, driveMag, [1, 3]);
angletUnique = unique(FFP_angle);
numAngle = length(angletUnique) ;


% close all
figure; plot(FFP_x, FFP_z); hold on; %  plot(FFP_x, z); 
xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2 FOV_z/2])
% k = 1.5e4;
% numSamplePerDrivePeriod = 1/f_drive*fs;
% idx = numSamplePerDrivePeriod*(k-1)+1:numSamplePerDrivePeriod*k;
% plot(FFP_x(idx), FFP_z(idx), 'linewidth', 5);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
load phantom1; % load phantom
tau = [2e-6 3e-6];
SPIOdistribution = zeros(512, 512, 2);
SPIOdistribution(1:40, (225:275), 1) = 1;
SPIOdistribution((1:40)+100, (225:275), 2) = 1;
% surf(SPIOdistribution(:,:, 1)); view(2); shading interp

image_FOV_x = 0.05;
image_FOV_z = 0.05;

dx = image_FOV_x/size(SPIOdistribution,1);  % distance between each pixel (m)
dz = image_FOV_z/size(SPIOdistribution,2); 

% calculate colinear and trasnverse PSF(s) for each unique angle
[colinearPSF, transversePSF, X, Z] = generatePSF([Gxx Gyy Gzz], 21, FOV_x, FOV_z, dx, dz, angletUnique, numAngle); 

   
% pad the empty parts of the image so that it has same size with PSF(s),
% then compute colinear and transverse images for each unique angle
for k=1:length(tau)
    [colinearIMG{k}, transverseIMG{k}] = generatePSFImages(SPIOdistribution(:,:,k), colinearPSF, transversePSF);
end


downsample = fs/fs_mpi;
for k=1:length(tau)
    % generate colinear, transverse, horizontal coil and vertical coil signals
    t_tau = (0:1/fs:tau(k)*15);
    r_t = 1/tau(k)*exp(-t_tau./tau(k));
    r_t = r_t/sum(r_t);
    [colinearSignal, transverseSignal, horizontalSignal{k}, verticalSignal{k}] = ...
        generateSignals(colinearIMG{k}, transverseIMG{k}, FFP_x, FFP_z, FFP_speed, FFP_angle, angletUnique, FOV_x, FOV_z, dx, dz, r_t);

    % downsample and high pass filter the received signals
    horizontalSignal_mpi{k} = horizontalSignal{k}(1:downsample:end);
    verticalSignal_mpi{k} = verticalSignal{k}(1:downsample:end);
end

horizontalSignal_mpi_mat = zeros(size(horizontalSignal_mpi{1}));
for k=1:length(tau)
    horizontalSignal_mpi_mat = horizontalSignal_mpi_mat + horizontalSignal_mpi{k};
end

interp_coeff = 2;
numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
numPeriods = length(horizontalSignal_mpi_mat)/numSamplePerDrivePeriod;
tau_est = zeros(1, numPeriods);
for k=1:numPeriods
    numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
    idx = (numSamplePerDrivePeriod*(k-1)+1:numSamplePerDrivePeriod*k);
    t_sig = idx/fs_mpi;
    
    numSamplePerDrivePeriod = 1/f_drive*(fs_mpi*interp_coeff);
    idx_interp = (numSamplePerDrivePeriod*(k-1)+1:numSamplePerDrivePeriod*k);
    t_interp = idx_interp/fs_mpi/interp_coeff;
   
    sig = interp1(t_sig, horizontalSignal_mpi_mat(idx),t_interp, 'spline');
    
    pos = sig(1:end/2); 
    neg = sig(end/2+1:end);
    S2=fft(neg);
    S1=fft(pos);
    sum_val = (S1+conj(S2));
    sub_val = (conj(S2)-S1);
    L = length(sum_val);
    f = (0:L-1)*(fs_mpi*interp_coeff)/L-(fs_mpi*interp_coeff)/2;
    f = fftshift(f);
    tau_est(k) = mean(real(sum_val(2:5)./(2*pi*1i*sub_val(2:5).*f(2:5))));
end

figure; scatter3(FFP_x(1:downsample*5:end-1), FFP_z(1:downsample*5:end-1), horizontalSignal_mpi_mat(1:5:end), 4,horizontalSignal_mpi_mat(1:5:end));  view(2); 
xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2 FOV_z/2])

figure; scatter3(FFP_x(1:2000:end-1), FFP_z(1:2000:end-1), tau_est, 4,tau_est);  view(2); 
xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2 FOV_z/2])


