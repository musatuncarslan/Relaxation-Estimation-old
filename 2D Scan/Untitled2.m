clear all
close all
clc

test = 1;

format long;

addpath('./parameterFiles')
addpath('./signalGenerationFiles')

% parameters
Physicsparams = setPhysicsParams(); % physics parameters
MPIparams = setMPIParams(Physicsparams); % MPI machine parameters
SPIOparams = setSPIOParams(Physicsparams); % SPIO parameters
Simparams = setSimulationParams(MPIparams, Physicsparams); % Simulation parameters


fileName = ['./signalSaveFolder/signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
% [fileObj, signal_size, chunk] = setFileParams(fileName, Physicsparams, MPIparams, SPIOparams, Simparams);
fileObj = matfile(fileName);
fileObj.Properties.Writable = true;

% create space for the simulation data
signal_size = Simparams.numSamplesPerIter/Simparams.downsample;
chunk = Simparams.numSamplesPerIter/Simparams.downsample;
numIters = Simparams.endIter-Simparams.startIter;
fileObj.horizontalSignal_mpi_mat(numIters,signal_size) = 0;
fileObj.FFP_x(numIters,signal_size) = 0;
fileObj.FFP_z(numIters,signal_size) = 0;
fileObj.FFP_speed(numIters,signal_size) = 0;
fileObj.FFP_angle(numIters,signal_size) = 0;


% start simulation
nout = 0;
count = 1;

for k = Simparams.startIter:Simparams.endIter
    tic
    
    t = ((k-1)*Simparams.numSamplesPerIter:k*Simparams.numSamplesPerIter)/Physicsparams.fs;
    
    FFPparams = generateFFP(t, MPIparams, Simparams);
   
    % calculate colinear and trasnverse PSF(s) for each unique angle
    
    [colinearPSF, transversePSF, X, Z] = generatePSF(MPIparams, SPIOparams, FFPparams);


    % pad the empty parts of the image so that it has same size with PSF(s),
    % then compute colinear and transverse images for each unique angle
    [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF, transversePSF);


    % generate MPI signals
    [signals] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams);


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
% 
% collinearSignal = horizontalSignal_mpi_mat.*cosd(FFP_angle) + verticalSignal_mpi_mat.*sind(FFP_angle);
% plot(FFP_z, collinearSignal./FFP_angle)
% 
% 
figure; scatter3(FFP_x, FFP_z, horizontalSignal_mpi_mat, 4, horizontalSignal_mpi_mat);  view(2); 
xlim([-MPIparams.FOV_x/2 MPIparams.FOV_x/2]); ylim([-MPIparams.FOV_z/2 MPIparams.FOV_z/2])