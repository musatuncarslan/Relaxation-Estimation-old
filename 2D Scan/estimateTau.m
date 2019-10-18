clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

filePath = 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder\';
% filePath = 'C:\Users\Orion\Google Drive\Phd\';
filename = [filePath, 'signal_18-October-2019_12-08-19'];
matObj = matfile(filename);

% get simulation parameters
MPIparams = matObj.MPIparams;
SPIOparams = matObj.SPIOparams;
Physicsparams = matObj.Physicsparams;
Simparams = matObj.Simparams;


f_d = MPIparams.f_drive; 
slewRate = MPIparams.slewRate; % T/s
Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
G=MPIparams.Gzz/Physicsparams.mu0; % gradient
MPIparams.driveMag=Hp/G; % extent of the drive field

R = slewRate/MPIparams.Gzz;

p = zeros(1, 21);
sign = zeros(1,ceil(length(p)/2));
for n=1:ceil(length(p)/2)
    sign(n) = (-1).^(n-1);
end

arg = 2*pi*f_d;

count = 1;
for k=1:length(p)
    if mod(k, 2) == 1
        p(k) = sign(count)*arg^k/factorial(k);
        count = count + 1;
    else
        p(k) = 0;
    end
end
p = MPIparams.driveMag*flip(p);
p_temp = p;
p_temp(end) = p_temp(end)+R;
roots_p = roots(p_temp);
real_idx = imag(roots_p) == 0;
real_roots = roots_p(real_idx);
del_t = real_roots(end) - 1/f_d/2;
 




Gxx = MPIparams.Gxx; % gradients (T/m)
Gyy = MPIparams.Gyy;
Gzz = MPIparams.Gzz;
Bp = MPIparams.Bp; % drive field (T)
f_drive = MPIparams.f_drive; % drive field frequency (Hz)
fs_mpi = MPIparams.fs; % sampling rate of the MPI machine (Hz)
FOV_x = MPIparams.FOV_x; % FOV in x-axis (m)
FOV_z = MPIparams.FOV_z; % FOV in z-axis (m)
time = MPIparams.time; % time on the line for slew rate calculation

diameter = SPIOparams.diameter; % nanoparticle diameter

fs = Physicsparams.fs; % sampling rate of the physics (Hz)

total_time = round(Simparams.numIters*Simparams.numSamplesPerIter/Physicsparams.fs, 5);

[numIters, numSamplesPerIter] = size(matObj.horizontalSignal_mpi_mat);
numIters = Simparams.numIters;
numPeriod = total_time*f_drive; % number of periods on a single line

interp_coeff = 600;
numPeriodsPerIter = numPeriod/numIters;
numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
pos_idx = (numSamplePerDrivePeriod/2:numSamplePerDrivePeriod:numSamplesPerIter)-50+1;
sig_contribution = 2:4;
count = 1;
tau_est_frequency = [];
tau_est_linear = [];
tau_lin_weighted = [];
FFP_x = [];
FFP_z = [];
for k=1:numIters
    tic
    numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
    idx = (numSamplePerDrivePeriod*numPeriodsPerIter*(k-1)+1:numSamplePerDrivePeriod*k*numPeriodsPerIter+1);
    t_sig = idx/fs_mpi;



    numSamplePerDrivePeriod = 1/f_drive*(fs_mpi*interp_coeff);
    idx_interp = (numPeriodsPerIter*numSamplePerDrivePeriod*(k-1)+1:numSamplePerDrivePeriod*numPeriodsPerIter*k+1);
    t_interp = idx_interp/fs_mpi/interp_coeff;

    sig = interp1(t_sig, matObj.horizontalSignal_mpi_mat(k, 1:numSamplesPerIter),t_interp, 'spline');
    for l=1:numPeriodsPerIter
%         periodIdx = (5*numSamplePerDrivePeriod*(l-1)+1:4.5*numSamplePerDrivePeriod*l+1);
%         partialSig = sig(periodIdx);
        partialSig = sig;

        L = (length(sig)-1-numSamplePerDrivePeriod/2);
        f = (0:L-1)*(fs_mpi*interp_coeff)/L-(fs_mpi*interp_coeff)/2;
        f = fftshift(f);
        
        pos = partialSig(1:end-1-numSamplePerDrivePeriod/2); 
        neg = partialSig(numSamplePerDrivePeriod/2+2:end);
        
        S2=fft(neg).*exp(-1i*2*pi*del_t.*f).*(Bp/Gzz*2*pi*f_d*cos(2*pi*f_d*del_t) - R)/(Bp/Gzz*2*pi*f_d + R);
        S1=fft(pos);
        
        
        sum_val = (S1+conj(S2));
        sub_val = (conj(S2)-S1);




        a = 1i*2*pi*f.*sub_val;
        b = sum_val;
        
%         A = toeplitz(real(ifft(a)));
%         A_inv = pinv(A);
%         b_t = transpose(real(ifft(b)));
%         tau = A_inv*b_t;
%         tau_t(count, l) = tau(1);
        
        a = a(sig_contribution);
        b = b(sig_contribution);

        tau_est_frequency(count, l) = mean(real(b./a));
        tau_est_linear(count, l) = real(transpose(pinv(a'*a)*a')*transpose(b));

        W = diag(abs(S1(sig_contribution)).^2);
        a_t = transpose(a);
        b_t = transpose(b);
        tau_lin_weighted(count, l) = real(inv(a_t'*W*a_t)*(W*a_t)'*b_t);
        
        
    end
    toc

    FFP_x(count, :) = matObj.FFP_x(k, pos_idx);
    if strcmp(MPIparams.ffp_type, 'linear_rastered')
        FFP_z(count, :) = matObj.FFP_z(k, pos_idx);
    else
        FFP_z(count, :) = mean(matObj.FFP_z(k, pos_idx))*ones(size(pos_idx));
    end

    count = count + 1;
end

x_axis = transpose(FFP_x(:, 1:end));
z_axis = transpose(FFP_z(:, 1:end));
tau_est_frequency = transpose(tau_est_frequency(:, 1:end))*1e6;
tau_est_linear = transpose(tau_est_linear(:, 1:end))*1e6;
tau_lin_weighted = transpose(tau_lin_weighted(:, 1:end))*1e6;
% tau_t = transpose(tau_t)*1e6;
% figure; scatter3(FFP_x(1:end), FFP_z(1:end), tau_est(1:end), 4,tau_est(1:end));  view(2); hold on;
% view(2); xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2, FOV_z/2]); xlabel('x-axis (m)'); ylabel('z-axis (m)'); zlabel('\tau'); colorbar;


SPIOdistribution = SPIOparams.SPIOdistribution;
% surf(SPIOdistribution(:,:, 1)); view(2); shading interp

tau = SPIOparams.tau*1e6;
tau_image = zeros(size(SPIOdistribution(:,:,1)));
for k=1:length(SPIOparams.diameter)
    tau_image = tau_image + SPIOdistribution(:, :, k)*tau(k);
end
img_size = size(tau_image); psf_size = Simparams.psf_size;
if img_size(1) < psf_size(1)
   tau_image = padarray(tau_image,[floor((psf_size(1)-img_size(1))/2), 0],0,'both');
   img_size = size(tau_image);
   tau_image = padarray(tau_image,[psf_size(1)-img_size(1), 0],0,'post');
end
if img_size(2) < psf_size(2)
   tau_image = padarray(tau_image,[0, floor((psf_size(2)-img_size(2))/2)],0,'both');
   img_size = size(tau_image);
   tau_image = padarray(tau_image,[0, psf_size(2)-img_size(2)],0,'post');
end

image_FOV_x = FOV_x;
image_FOV_z = FOV_z;

lin_z_axis = linspace(-image_FOV_z/2, image_FOV_z/2, size(tau_image, 1));
z_vals = lin_z_axis(find(tau_image(:, 256) == 2));


fig = figure('Position', [560 240 800 600]);
plot(z_axis(:), tau_est_frequency(:), 'linewidth', 2, 'marker', 'o'); hold on; 
plot(z_axis(:), tau_est_linear(:), 'linewidth', 2, 'marker', 'd'); hold on; 
plot(z_axis(:), tau_lin_weighted(:), 'linewidth', 2, 'marker', '*'); hold on; 
% plot(z_axis(:), tau_t(:), 'linewidth', 2, 'marker', '*'); hold on; 
plot(lin_z_axis, tau_image(:, 256), 'linewidth', 2)
xlabel('z-axis (m)'); ylabel('\tau (\mu s)'); axis tight; 
legend('Frequency Estimated \tau', 'Linear Estimated \tau', 'Weighted Linear Estimated \tau', 'Original \tau') % ,'t Domain Deconvolution'

if strcmp(MPIparams.ffp_type, 'linear_rastered')
    type_str = 'Rastered Linear';
    slew_rate_str = ['R_s = ' num2str(MPIparams.slewRate ) ' T/s'];
else
    type_str = 'Fixed FFP';
    slew_rate_str = ['R_S = n/a'];
end
gradent_str = ['G_x = ' num2str(Gxx), ', G_y = ' num2str(Gyy), ', G_z = ' num2str(Gzz) ' T/m'];
diameter_str = ['Diameter = ' num2str(diameter) ' nm'];
Bp_str = ['Bp = ' num2str(Bp*1e3) ' mT'];
f_drive_str = ['F_{drive} = ' num2str(f_drive*1e-3) ' kHz'];
fs_str = ['F_s = ' num2str(fs*1e-6) ' MHz'];
fs_mpi_str = ['F_s_{mpi} = ' num2str(fs_mpi*1e-6) ' MHz'];

str = {type_str, '% Parameters', gradent_str, diameter_str, Bp_str, f_drive_str, ...
    fs_str, fs_mpi_str, slew_rate_str};

legend_pos = fig.Children(1).Position;
pos = [legend_pos(1), 0.2, 0.4 0.3];
annotation('textbox', pos, 'String',str,'FitBoxToText','on');
% ylim([1.85 2]); xlim([-0.03 0.03])
