clear all
close all
clc

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

filePath = 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder\';
% filePath = 'C:\Users\Orion\Google Drive\Phd\';
filename = [filePath, 'signal_29-October-2019_16-51-01'];
matObj = matfile(filename);

% get simulation parameters
MPIparams = matObj.MPIparams;
SPIOparams = matObj.SPIOparams;
Physicsparams = matObj.Physicsparams;
Simparams = matObj.Simparams;

[del_t, a_t] = findDistortion(MPIparams); 
del_t = del_t - 1/2/MPIparams.f_drive;
a_t = -a_t;

Gzz = MPIparams.Gzz;
Bp = MPIparams.Bp; % drive field (T)
f_drive = MPIparams.f_drive; % drive field frequency (Hz)
fs_mpi = MPIparams.fs; % sampling rate of the MPI machine (Hz)
FOV_z = MPIparams.FOV_z; % FOV in z-axis (m)

traversedFOVz = MPIparams.traversedFOVz;
robotSpeed = MPIparams.slewRate/MPIparams.Gzz; % robot arm movement speed (m/s)
total_time = (traversedFOVz(2)-traversedFOVz(1))/robotSpeed;


numIters = Simparams.numIters;
numPeriod = Simparams.numIters*Simparams.numSamplesPerIter*MPIparams.f_drive/Physicsparams.fs; % number of periods on a single line

interp_coeff = 100;
numPeriodsPerIter = numPeriod/numIters;
numSamplesPerIter = numPeriodsPerIter*(1/f_drive*fs_mpi) + 2;
numSample = (1/f_drive*fs_mpi) + 2; % + 2 is for interpolation reasons
numSampleInterpolated = numSample*interp_coeff;

sig_contribution = 2:4;
count = 1;
tau_est_frequency = [];
tau_est_linear = [];
tau_lin_weighted = [];
for k=1:numIters
    tic
    idx = 1:numSamplesPerIter;
    t_sig = idx/fs_mpi;

    idx_interp = (1:numSampleInterpolated);
    t_interp = idx_interp/fs_mpi/interp_coeff;

    sig = interp1(t_sig, matObj.horizontalSignal_mpi_mat(k, 1:numSamplesPerIter),t_interp, 'spline');
    for l=1:numPeriodsPerIter
        partialSig = sig;

        pos = partialSig(1:MPIparams.fs/MPIparams.f_drive/2*interp_coeff); 
        neg = partialSig(MPIparams.fs/MPIparams.f_drive/2*interp_coeff+2*interp_coeff:2*MPIparams.fs/MPIparams.f_drive/2*interp_coeff+2*interp_coeff-1);

        L = length(neg);
        f = (0:L-1)*(fs_mpi*interp_coeff)/L-(fs_mpi*interp_coeff)/2;
        f = fftshift(f);

        S1 = fft(pos);
        S2=fft(neg).*exp(1i*2*pi*(round(del_t-1/fs_mpi/interp_coeff/5, 9)).*f);
%         shift_func = real(ifft(exp(1i*2*pi*(round(del_t, 9)).*f)));
%         shift_func(shift_func ~= 1) = 0;
%         s2 = conv(neg, shift_func);
%         s2 = s2(1:MPIparams.fs/MPIparams.f_drive/2*interp_coeff).*a_t;      
%         S2 = fft(s2);



        sum_val = (S1+conj(S2));
        sub_val = (conj(S2)-S1);

        a = 1i*2*pi*f.*sub_val;
        b = sum_val;

        plot(real(b./a))

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

    count = count + 1;
end

tau_est_frequency = transpose(tau_est_frequency(:, 1:end))*1e6;
tau_est_linear = transpose(tau_est_linear(:, 1:end))*1e6;
tau_lin_weighted = transpose(tau_lin_weighted(:, 1:end))*1e6;

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

image_FOV_z = FOV_z;

lin_z_axis = linspace(-image_FOV_z/2, image_FOV_z/2, size(tau_image, 1));
z_vals = lin_z_axis(find(tau_image(:, 512) == 2));



plot(z_vals, tau_est_frequency(:), 'linewidth', 2, 'marker', 'o', 'color', [0 0.4470 0.7410]); hold on; 
plot(z_vals, tau_est_linear(:), 'linewidth', 2, 'marker', 'd', 'color', [0.8500 0.3250 0.0980]); hold on; 
plot(z_vals, tau_lin_weighted(:), 'linewidth', 2, 'marker', '*', 'color', [0.9290 0.6940 0.1250]); hold on; 
% plot(z_axis(:), tau_t(:), 'linewidth', 2, 'marker', '*'); hold on; 
if (plotif == 1)
    plot(lin_z_axis, tau_image(:, 512), 'linewidth', 2, 'color', [0.4940 0.1840 0.5560])
end
xlabel('z-axis (m)'); ylabel('\tau (\mu s)'); axis tight; 