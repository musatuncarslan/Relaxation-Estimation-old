clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 20); % MPI machine parameters

fig = figure('Position', [560 240 800 600]); hold on;
pos = 512;
for pos_idx=1:length(pos)
    SPIOparams = setSPIOParams(Physicsparams, pos(pos_idx), 2e-6); % SPIO parameters
    Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

    filePath = 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder\';
    fileName = ['signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
    fullFileName = [filePath fileName];
    % [fileObj, signal_size, chunk] = setFileParams(fileName, Physicsparams, MPIparams, SPIOparams, Simparams);
    fileObj = matfile(fullFileName);
    fileObj.Properties.Writable = true;

    % create space for the simulation data
    signal_size = Simparams.numSamplesPerIter/Simparams.downsample+1;
    chunk = Simparams.numSamplesPerIter/Simparams.downsample;
    numIters = Simparams.numIters;
    fileObj.horizontalSignal_mpi_mat(numIters,signal_size) = 0;
    fileObj.FFP_x(numIters,signal_size) = 0;
    fileObj.FFP_z(numIters,signal_size) = 0;
    fileObj.FFP_speed(numIters,signal_size)  = 0;
    fileObj.FFP_angle(numIters,signal_size) = 0;


    % start simulation
    nout = 0;
    count = 1;
    for k = Simparams.startIter:Simparams.endIter
        tic

        t = ((k-1)*Simparams.numSamplesPerIter+1:(k*Simparams.numSamplesPerIter+Simparams.downsample+1))/Physicsparams.fs;
        FFPparams = generateFFP(t, MPIparams, Simparams);

        % calculate colinear and trasnverse PSF(s) for each unique angle

        [colinearPSF, transversePSF, X, Z] = generatePSF(MPIparams, SPIOparams, FFPparams);


        % pad the empty parts of the image so that it has same size with PSF(s),
        % then compute colinear and transverse images for each unique angle
        [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF, transversePSF);


        % generate MPI signals
        [signals_sep] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams);

        signals = struct;
        signals.horizontalSignal = zeros(1, size(signals_sep, 2));
        signals.verticalSignal = zeros(1, size(signals_sep, 2));
        for l=1:length(SPIOparams.diameter)
            signals.horizontalSignal = signals.horizontalSignal + signals_sep.horizontalSignal(l, :);
            signals.verticalSignal = signals.verticalSignal + signals_sep.verticalSignal(l, :);
        end

        % save signals to the file
        fprintf('Writing %d of %d, ',nout,signal_size*Simparams.numIters);
        chunkSize = min(chunk,signal_size*Simparams.numIters-nout);
        fileObj.horizontalSignal_mpi_mat(count,1:signal_size) = signals.horizontalSignal;
        fileObj.verticalSignal_mpi_mat(count,1:signal_size) = signals.verticalSignal;
        fileObj.FFP_x(count,1:signal_size) = FFPparams.FFP_x(1:Simparams.downsample:end-1);
        fileObj.FFP_z(count,1:signal_size) = FFPparams.FFP_z(1:Simparams.downsample:end-1);
        fileObj.FFP_speed(count,1:signal_size) = FFPparams.FFP_speed(1:Simparams.downsample:end);
        fileObj.FFP_angle(count,1:signal_size) = FFPparams.FFP_angle(1:Simparams.downsample:end);
        nout = nout + chunkSize;

        endTime(count) = toc;
        fprintf('Iter: %i of %i, Time: %d\n', count, Simparams.numIters, endTime(count));
        count = count + 1;
    end

    totalTime = sum(endTime)

    Simparams.psf_size = size(colinearPSF{1, 1});
    % save simulation parameters
    fileObj.MPIparams = MPIparams;
    fileObj.SPIOparams = SPIOparams;
    fileObj.Physicsparams = Physicsparams;
    fileObj.Simparams = Simparams;


    FFP_x = transpose(fileObj.FFP_x); FFP_x = FFP_x(:);
    FFP_z = transpose(fileObj.FFP_z); FFP_z = FFP_z(:);
    FFP_speed = transpose(fileObj.FFP_speed); FFP_speed = FFP_speed(:);
    FFP_angle = transpose(fileObj.FFP_angle); FFP_angle = FFP_angle(:);
    horizontalSignal_mpi_mat = transpose(fileObj.horizontalSignal_mpi_mat); horizontalSignal_mpi_mat = horizontalSignal_mpi_mat(:);
    verticalSignal_mpi_mat = transpose(fileObj.verticalSignal_mpi_mat); verticalSignal_mpi_mat = verticalSignal_mpi_mat(:);

        
    if ((pos_idx == 1) || pos_idx == length(pos))
        plotif = 1;
    else
        plotif = 0;
    end
    estimateTauFunc(fileName, plotif)

end
plot([min(FFPparams.FFP_z), max(FFPparams.FFP_z)], ones(1, 2), 'linewidth', 2, 'color', [0.4660 0.6740 0.1880]);
if strcmp(MPIparams.ffp_type, 'linear_rastered')
    type_str = 'Rastered Linear';
    slew_rate_str = ['R_s = ' num2str(MPIparams.slewRate) ' T/s'];
else
    type_str = 'Fixed FFP';
    slew_rate_str = 'R_S = n/a';
end
gradent_str = ['G_x = ' num2str(MPIparams.Gxx), ', G_y = ' num2str(MPIparams.Gyy), ', G_z = ' num2str(MPIparams.Gzz) ' T/m'];
diameter_str = ['Diameter = ' num2str(SPIOparams.diameter) ' nm'];
Bp_str = ['Bp = ' num2str(MPIparams.Bp*1e3) ' mT'];
f_drive_str = ['F_{drive} = ' num2str(MPIparams.f_drive*1e-3) ' kHz'];
fs_str = ['F_s = ' num2str(Physicsparams.fs*1e-6) ' MHz'];
fs_mpi_str = ['F_s_{mpi} = ' num2str(MPIparams.fs*1e-6) ' MHz'];

str = {type_str, '% Parameters', gradent_str, diameter_str, Bp_str, f_drive_str, ...
    fs_str, fs_mpi_str, slew_rate_str};

xlabel('z-axis (m)'); ylabel('\tau (\mu s)'); axis tight; 
legend_pos = fig.Children(1).Position;
pos = [legend_pos(1), 0.2, 0.4 0.3];
annotation('textbox', pos, 'String',str,'FitBoxToText','on');
h = findobj(gca,'Type','line');
legend([h(5) h(4) h(3) h(2) h(1)], {'Frequency Estimated \tau', 'Linear Estimated \tau', ...
    'Weighted Linear Estimated \tau', 'Original \tau', 'Extent of FFP'}) % ,'t Domain Deconvolution'
