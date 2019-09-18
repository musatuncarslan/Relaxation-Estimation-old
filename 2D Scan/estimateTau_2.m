aaaclear all
close all
clc

format long;

filenames = dir('*.mat');
filenames = filenames(end-17:end-2);

for m=1:length(filenames)

    matObj = matfile(filenames(m).name);

    % get simulation parameters
    Gxx = matObj.params(1, 1); % gradients (T/m)
    Gyy = matObj.params(1, 2);
    Gzz = matObj.params(1, 3);
    diameter = matObj.params(1, 4); % nanoparticle diameter
    Bp = matObj.params(1, 5); % drive field (T)
    f_drive = matObj.params(1, 6); % drive field frequency (Hz)
    fs = matObj.params(1, 7); % sampling rate of the physics (Hz)
    fs_mpi = matObj.params(1, 8); % sampling rate of the MPI machine (Hz)
    FOV_x = matObj.params(1, 9); % FOV in x-axis (m)
    FOV_z = matObj.params(1, 10); % FOV in z-axis (m)
    time = matObj.params(1, 11); % time robot arm moves (s), this has double 
                                 % meaning, for zig-zag movement, its the total
                                 % time it takes for robot arm to traverse 
                                 % whole FOV in 2D, for 1D movement, its the
                                 % time robot arm needs to traverse 1 line.



    [numIters, numSamplePerIter] = size(matObj.horizontalSignal_mpi_mat);
    numPeriod = time*f_drive; % number of periods on a single line

    interp_coeff = 1;
    numPeriodsPerIter = numPeriod/numIters;
    numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
    pos_idx = (numSamplePerDrivePeriod/2:numSamplePerDrivePeriod:numSamplePerIter);
    sig_contribution = 2:5;
    count = 1;
    tau_est_frequency = [];
    tau_est_linear = [];
    tau_lin_weighted = [];
    FFP_x = [];
    FFP_z = [];
    for k=1:numIters
        tic
        numSamplePerDrivePeriod = 1/f_drive*fs_mpi;
        idx = (numSamplePerDrivePeriod*numPeriodsPerIter*(k-1)+1:numSamplePerDrivePeriod*k*numPeriodsPerIter);
        t_sig = idx/fs_mpi;



        numSamplePerDrivePeriod = 1/f_drive*(fs_mpi*interp_coeff);
        idx_interp = (numPeriodsPerIter*numSamplePerDrivePeriod*(k-1)+1:numSamplePerDrivePeriod*numPeriodsPerIter*k);
        t_interp = idx_interp/fs_mpi/interp_coeff;

        sig = interp1(t_sig, matObj.horizontalSignal_mpi_mat(k, 1:numSamplePerIter),t_interp, 'spline');
        for l=1:numPeriodsPerIter
            periodIdx = (numSamplePerDrivePeriod*(l-1)+1:numSamplePerDrivePeriod*l);
            partialSig = sig(periodIdx);

            pos = partialSig(1:end/2); 
            neg = partialSig(end/2+1:end);
            S2=fft(neg);
            S1=fft(pos);
            sum_val = (S1+conj(S2));
            sub_val = (conj(S2)-S1);


            L = length(sum_val);
            f = (0:L-1)*(fs_mpi*interp_coeff)/L-(fs_mpi*interp_coeff)/2;
            f = fftshift(f);

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
        if strcmp(matObj.type, 'linear_rastered')
            FFP_z(count, :) = matObj.FFP_z(k, pos_idx);
        else
            FFP_z(count, :) = mean(matObj.FFP_z)*ones(size(pos_idx));
        end

        count = count + 1;
    end

    x_axis = transpose(FFP_x(:, 2:end));
    z_axis = transpose(FFP_z(:, 2:end));
    tau_est_frequency = transpose(tau_est_frequency(:, 2:end))*1e6;
    tau_est_linear = transpose(tau_est_linear(:, 2:end))*1e6;
    tau_lin_weighted = transpose(tau_lin_weighted(:, 2:end))*1e6;
    % figure; scatter3(FFP_x(1:end), FFP_z(1:end), tau_est(1:end), 4,tau_est(1:end));  view(2); hold on;
    % view(2); xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2, FOV_z/2]); xlabel('x-axis (m)'); ylabel('z-axis (m)'); zlabel('\tau'); colorbar;


    SPIOdistribution = zeros(512, 512, 2);
    SPIOdistribution(1:40, (225:275), 1) = 1;
    SPIOdistribution((1:40)+100, (225:275), 2) = 1;
    % surf(SPIOdistribution(:,:, 1)); view(2); shading interp

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

%     f = figure('Position', [560 240 800 600]);
%     plot(z_axis(:), tau_est_frequency(:), 'linewidth', 2, 'marker', 'o'); hold on; 
%     plot(z_axis(:), tau_est_linear(:), 'linewidth', 2, 'marker', 'd'); hold on; 
%     plot(z_axis(:), tau_lin_weighted(:), 'linewidth', 2, 'marker', '*'); hold on; 
%     plot(lin_z_axis, tau_image(:, 256), 'linewidth', 2)
%      xlabel('z-axis (m)'); ylabel('\tau (\mu s)'); axis tight; legend('Frequency Estimated \tau', 'Linear Estimated \tau', 'Weighted Linear Estimated \tau', 'Original \tau')
% 
%     if strcmp(matObj.type, 'linear_rastered')
%         type_str = 'Linear Rastered';
%     else
%         type_str = 'Fixed';
%     end
%     gradent_str = ['G_x = ' num2str(Gxx), ', G_y = ' num2str(Gyy), ', G_z = ' num2str(Gzz) ' T/m'];
%     diameter_str = ['Diameter = ' num2str(diameter) ' nm'];
%     Bp_str = ['Bp = ' num2str(Bp*1e3) ' mT'];
%     f_drive_str = ['F_{drive} = ' num2str(f_drive*1e-3) ' kHz'];
%     fs_str = ['F_s = ' num2str(fs*1e-6) ' MHz'];
%     fs_mpi_str = ['F_s_{mpi} = ' num2str(fs_mpi*1e-6) ' MHz'];
%     slew_rate_str = ['R_s = ' num2str(abs(FOV_z/time*Gzz)) ' T/s'];
% 
%     str = {type_str, '% Parameters', gradent_str, diameter_str, Bp_str, f_drive_str, ...
%         fs_str, fs_mpi_str, slew_rate_str};
% 
%     legend_pos = f.Children(1).Position;
%     pos = [legend_pos(1), 0.2, 0.4 0.3];
%     annotation('textbox', pos, 'String',str,'FitBoxToText','on');

end

figure
[B,I] = sort(tau_estimations.slew_rate(2:end-1));
plot(B, tau_estimations.tau_freq(I+1), '-o', 'linewidth', 2); hold on;
plot(B, tau_estimations.tau_linear(I+1), '-d', 'linewidth', 2)
plot(B, tau_estimations.tau_weighted_linear(I+1), '-*', 'linewidth', 2)
plot(B, ones(1, length(tau_estimations.slew_rate(2:end-1)))*tau_estimations.tau_weighted_linear(end), '-', 'linewidth', 2)
xlabel('Slew Rate (T/s)'); ylabel('\tau (s)')
legend('Frequency Estimation', 'Least Squares Solution', 'Weighted Least Squares', 'Fixed Robot Arm (Freq. Est.)')