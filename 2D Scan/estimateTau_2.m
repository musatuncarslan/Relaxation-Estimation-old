clear all
close all
clc

format long;

filePath = 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder';
filenames = dir(filePath);
filenames = filenames(3:end);

for m=1:length(filenames)

    matObj = matfile([filePath, '\', filenames(m).name]);

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

    interp_coeff = 1;
    numPeriodsPerIter = numPeriod/numIters;
    numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
    pos_idx = (numSamplePerDrivePeriod/2:numSamplePerDrivePeriod:numSamplesPerIter)-50;
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
            periodIdx = (numSamplePerDrivePeriod*(l-1)+1:numSamplePerDrivePeriod*l+1);
            partialSig = sig(periodIdx);

            L = length(partialSig)/2;
            f = (0:L-1)*(fs_mpi*interp_coeff)/L-(fs_mpi*interp_coeff)/2;
            f = fftshift(f);


            pos = partialSig(1:floor(end/2)); 
            neg = partialSig(floor(end/2)+2:end);
            S2=fft(neg);
            S1=fft(pos);
            sum_val = (S1+conj(S2));
            sub_val = (conj(S2)-S1);




            a = 1i*2*pi*f(sig_contribution).*sub_val(sig_contribution);
            b = sum_val(sig_contribution);

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
    x_axis = transpose(FFP_x(:, :));
    z_axis = transpose(FFP_z(:, :));
    tau_est_frequency = transpose(tau_est_frequency(:, :))*1e6;
    tau_est_linear = transpose(tau_est_linear(:, :))*1e6;
    tau_lin_weighted = transpose(tau_lin_weighted(:, :))*1e6;
    % figure; scatter3(FFP_x(1:end), FFP_z(1:end), tau_est(1:end), 4,tau_est(1:end));  view(2); hold on;
    % view(2); xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2, FOV_z/2]); xlabel('x-axis (m)'); ylabel('z-axis (m)'); zlabel('\tau'); colorbar;


    SPIOdistribution = SPIOparams.SPIOdistribution;

    tau = [2 3];
    tau_image = SPIOdistribution(:, :, 1)*tau(1);

    img_size = size(tau_image); psf_size = [615, 513];
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

    [~, min_idx_1] = min(abs(z_axis(:)-z_vals(1)));
    [~, min_idx_2] = min(abs(z_axis(:)-z_vals(end)));

    tau_estimations.slew_rate(1, m) = abs(FOV_z/time*Gzz);
    tau_estimations.tau_freq(1, m) = mean(tau_est_frequency(min_idx_1:min_idx_2));
    tau_estimations.tau_linear(1, m) = mean(tau_est_linear(min_idx_1:min_idx_2));
    tau_estimations.tau_weighted_linear(1, m) = mean(tau_lin_weighted(min_idx_1:min_idx_2));



end

figure
[B,I] = sort(tau_estimations.slew_rate(1:end));
plot(B, tau_estimations.tau_freq(I), '-o', 'linewidth', 2); hold on;
plot(B, tau_estimations.tau_linear(I), '-d', 'linewidth', 2)
plot(B, tau_estimations.tau_weighted_linear(I), '-*', 'linewidth', 2)
plot(B, ones(1, length(B))*SPIOparams.tau(1)*1e6, 'linewidt', 2)
xlabel('Slew Rate (T/s)'); ylabel('\tau (\mu s)')
legend('Frequency Estimation', 'Least Squares Solution', 'Weighted Least Squares', 'Actual \tau Value')