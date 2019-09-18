clear all
close all
clc

format long;

robot_movement_type = 'fixed';
MPIparams = struct;
SPIOparams = struct;


% gradients (T/m) (current scanner)
Gxx = 4.8;
Gyy = 2.4;
Gzz = 2.4;
diameter = 25; % (nm)

Bp = 15e-3; % Drive field (T)
f_drive = 10e3; % drive field frequency
mu0=1.256637*10^-6; % permaebility of vacuum
Hp=Bp/mu0; % magnetization moment
G=Gzz/mu0; % gradient
driveMag=Hp/G; % extent of the drive field

fs = 20e6; % sample frequency of the physical world (Hz)
fs_mpi = 2e6; % sample frequency of the MPI system (Hz)
downsample = fs/fs_mpi; % downsample ratio
FOV_z = 0.06; % FOV in z-axis (meters) (bore axis)
FOV_x = 0.05; % FOV in x-axis (meters)
time = 2; % time (seconds)

numLines = 9; % number of lines in x-axis to traverse

noLine = 5;

numPeriod = time*f_drive; % number of periods on a single line
numSamplePerPeriod = fs/f_drive; % sample per drive field period in Fs
a = gcd(fs*time, numPeriod); % find the greatest common divisor of total number of samples and total number of periods per line
div = divisors(a); % find the divisors of gcd (these numbers divide both total number of samples and number of periods per line)
numIters = div(end-3); % use the middle number of divisors as a rule of thumb for number of iterations to solve the whole problem

numPeriodsPerIter = numPeriod/numIters;
numSamplesPerIter = fs*time/numIters;


tau = [2e-6];
SPIOdistribution = zeros(512, 512, 2);
SPIOdistribution(1:40, (225:275), 1) = 1;
SPIOdistribution((1:40)+100, (225:275), 2) = 1;
image_FOV_x = 0.05;
image_FOV_z = 0.05;
dx = image_FOV_x/size(SPIOdistribution,1);  % distance between each pixel (m)
dz = image_FOV_z/size(SPIOdistribution,2); 

for k=1:length(tau)
    [row,col] = find(SPIOdistribution(:, :, k) == 1);
    zloc =   unique(row)*dz - image_FOV_z/2;
    xloc = unique(col)*dx - image_FOV_x/2;
    SPIOlocation(k, :) = [xloc(ceil(end/2)) zloc(ceil(end/2))];
end

fileName = ['signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
matObj = matfile(fileName);
matObj.Properties.Writable = true;

% save simulation parameters
matObj.type = robot_movement_type;
matObj.params(1, 1:3) = [Gxx, Gyy, Gzz];
matObj.params(1, 4) = diameter;
matObj.params(1, 5) = Bp;
matObj.params(1, 6) = f_drive;
matObj.params(1, 7:8) = [fs, fs_mpi];
matObj.params(1, 9:11) = [FOV_x, FOV_z, time];

% create space for the simulation data
signal_size = numSamplesPerIter/downsample;
chunk = numSamplesPerIter/downsample;
matObj.horizontalSignal_mpi_mat(numIters,signal_size) = 0;
matObj.FFP_x(numIters,signal_size) = 0;
matObj.FFP_z(numIters,signal_size) = 0;
matObj.FFP_speed(numIters,signal_size) = 0;
matObj.FFP_angle(numIters,signal_size) = 0;

% start simulation
nout = 0;
x = linspace(-FOV_x/2, FOV_x/2, numLines); % robot arm movement in x direction w.r.t. time
for k = 1:numIters
    tic
    
    t = ((k-1)*numSamplesPerIter:k*numPeriodsPerIter*numSamplePerPeriod)/fs;
    
    if strcmp(robot_movement_type, 'linear_rastered')
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
    else
        z = repmat(SPIOlocation(1, 2), [1, numPeriodsPerIter*numSamplePerPeriod+1]);
    end
    drive = driveMag*cos(2*pi*f_drive*t); % drive field movement

    FFP_x = repmat(x(noLine), [1, numPeriodsPerIter*numSamplePerPeriod+1]); % movement of FFP in x direction, for linear x is constant throughout the line
    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)

    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFP_speed = sqrt(zDifference.^2 + xDifference.^2);
    FFP_angle = wrapTo360(round(atan2d(xDifference, zDifference), 1));    

%     [FFP_x, FFP_z, FFP_speed, FFP_angle, x, z] = FFPtrajectory(FOV_x, FOV_z, fs, f_drive, time, driveMag, 'triangular');
    angletUnique = unique(FFP_angle);
    numAngle = length(angletUnique) ;

   
    % calculate colinear and trasnverse PSF(s) for each unique angle
    [colinearPSF, transversePSF, X, Z] = generatePSF([Gxx Gyy Gzz], 21, FOV_x, FOV_z, dx, dz, angletUnique, numAngle); 


    % pad the empty parts of the image so that it has same size with PSF(s),
    % then compute colinear and transverse images for each unique angle
    for l=1:length(tau)
        [colinearIMG{l}, transverseIMG{l}] = generatePSFImages(SPIOdistribution(:,:,l), colinearPSF, transversePSF);
    end


    
    for l=1:length(tau)
        % generate colinear, transverse, horizontal coil and vertical coil signals
        t_tau = (0:1/fs:tau(l)*15);
        r_t = 1/tau(l)*exp(-t_tau./tau(l));
        r_t = r_t/sum(r_t);
        [colinearSignal, transverseSignal, horizontalSignal{l}, verticalSignal{l}] = ...
            generateSignals(colinearIMG{l}, transverseIMG{l}, FFP_x, FFP_z, FFP_speed, FFP_angle, angletUnique, FOV_x, FOV_z, dx, dz, r_t);

        % downsample and high pass filter the received signals
        horizontalSignal_mpi{l} = horizontalSignal{l}(1:downsample:end);
        verticalSignal_mpi{l} = verticalSignal{l}(1:downsample:end);
    end

    horizontalSignal_mpi_mat = zeros(size(horizontalSignal_mpi{1}));
    for l=1:length(tau)
        horizontalSignal_mpi_mat = horizontalSignal_mpi_mat + horizontalSignal_mpi{l};
    end
  
    fprintf('Writing %d of %d, ',nout,signal_size*numIters);
    chunkSize = min(chunk,signal_size*numIters-nout);
    matObj.horizontalSignal_mpi_mat(k,1:signal_size) = horizontalSignal_mpi_mat;
    matObj.FFP_x(k,1:signal_size) = FFP_x(1:downsample:end-1);
    matObj.FFP_z(k,1:signal_size) = FFP_z(1:downsample:end-1);
    matObj.FFP_speed(k,1:signal_size) = FFP_speed(1:downsample:end);
    matObj.FFP_angle(k,1:signal_size) = FFP_angle(1:downsample:end);
    nout = nout + chunkSize;
    
    
    
    endTime(k) = toc;
    fprintf('Iter: %i of %i, Time: %d\n', k, numIters, endTime(k));
end

totalTime = sum(endTime)

FFP_x = transpose(matObj.FFP_x); FFP_x = FFP_x(:);
FFP_z = transpose(matObj.FFP_z); FFP_z = FFP_z(:);
horizontalSignal_mpi_mat = transpose(matObj.horizontalSignal_mpi_mat); horizontalSignal_mpi_mat = horizontalSignal_mpi_mat(:);

figure; scatter3(FFP_x, FFP_z, horizontalSignal_mpi_mat, 4, horizontalSignal_mpi_mat);  view(2); 
xlim([-FOV_x/2 FOV_x/2]); ylim([-FOV_z/2 FOV_z/2])