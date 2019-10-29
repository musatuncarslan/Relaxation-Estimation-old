clear all
close all
clc

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams, 0); % MPI machine parameters

SPIOparams = setSPIOParams(Physicsparams, 512, 2e-6); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters

filePath = 'C:\Users\Orion\Google Drive\Phd\signalSaveFolder\';
fileName = ['signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
fullFileName = [filePath fileName];
[fileObj, signal_size, ~] = setFileParams(fullFileName, Physicsparams, MPIparams, SPIOparams, Simparams);

chunk = signal_size;
count = 1;
for k=Simparams.startIter:Simparams.endIter
    tic
    t = ((k-1)*Simparams.numSamplesPerIter+1:(k*Simparams.numSamplesPerIter+2*Simparams.downsample+1))/Physicsparams.fs;
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
    fileObj.horizontalSignal_mpi_mat(count,1:signal_size) = signals.horizontalSignal;
    fileObj.verticalSignal_mpi_mat(count,1:signal_size) = signals.verticalSignal;
    fileObj.FFP_x(count,1:signal_size) = FFPparams.FFP_x(1:Simparams.downsample:end-1);
    fileObj.FFP_z(count,1:signal_size) = FFPparams.FFP_z(1:Simparams.downsample:end-1);
    fileObj.FFP_speed(count,1:signal_size) = FFPparams.FFP_speed(1:Simparams.downsample:end);
    fileObj.FFP_angle(count,1:signal_size) = FFPparams.FFP_angle(1:Simparams.downsample:end);
    
    endTime(count) = toc;
    fprintf('%d of %d samples done in %2.4f sec.\n', chunk, signal_size*Simparams.numIters, endTime(count));
    chunk = chunk + signal_size;
    count = count + 1;
    
    
    
end
