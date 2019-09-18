clear all
close all
clc

format long;

MPIparams = struct;
SPIOparams = struct;

robot_movement_type = 'triangular_rastered';
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

fs = 0.5e6; % sample frequency of the physical world (Hz)
fs_mpi = 0.5e6; % sample frequency of the MPI system (Hz)
downsample = fs/fs_mpi; % downsample ratio
FOV_z = 0.05; % FOV in z-axis (meters) (bore axis)
FOV_x = 0.05; % FOV in x-axis (meters)
time = 0.5; % time (seconds)
traversedFOVz = [-0.025 0.025]; % traversed fov in the simulation in z-axis (m)

robotSpeed = FOV_z/time; % robot arm movement speed (m/s)
bFOVz = -FOV_z/2; % beginning point of FOV in z-axis (m)
t1 = round((traversedFOVz(1)-bFOVz)/robotSpeed, 2);
t2 = round((traversedFOVz(2)-bFOVz)/robotSpeed, 2);
total_time = round(t2-t1, 2);

numLines = 9; % number of lines in x-axis to traverse
noLine = 5;

numPeriod = total_time*f_drive; % number of periods on the portion of a single line
numSamplePerPeriod = fs/f_drive; % sample per drive field period in Fs
a = gcd(fs*total_time, numPeriod); % find the greatest common divisor of total number of samples and total number of periods per line
div = divisors(a); % find the divisors of gcd (these numbers divide both total number of samples and number of periods per line)
numIters = div(end-6); % use the middle number of divisors as a rule of thumb for number of iterations to solve the whole problem

numPeriodsPerIter = numPeriod/numIters;
numSamplesPerIter = fs*total_time/numIters;

startIter = (t1 /numSamplesPerIter)*fs+1;
endIter = startIter+numIters-1;

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

fileName = 'trajectory_test_file'; % ['signal_' datestr(now,'dd-mmmm-yyyy_HH-MM-SS') '.mat'];
matObj = matfile(fileName);
matObj.Properties.Writable = true;

% save simulation parameters
matObj.type = robot_movement_type;
matObj.params(1, 1:3) = [Gxx, Gyy, Gzz];
matObj.params(1, 4) = diameter;
matObj.params(1, 5) = Bp;
matObj.params(1, 6) = f_drive;
matObj.params(1, 7:8) = [fs, fs_mpi];
matObj.params(1, 9:13) = [FOV_x, FOV_z, t1, t2, time];

% create space for the simulation data
signal_size = numSamplesPerIter/downsample;
chunk = numSamplesPerIter/downsample;
matObj.horizontalSignal_mpi_mat(numIters,signal_size) = 0;
matObj.FFP_x(numIters,signal_size) = 0;
matObj.FFP_z(numIters,signal_size) = 0;
matObj.z(numIters, signal_size) = 0;
matObj.FFP_speed(numIters,signal_size) = 0;
matObj.FFP_angle(numIters,signal_size) = 0;

% start simulation
nout = 0;
count = 1;

for k = startIter:endIter
    tic
    
    t = ((k-1)*numSamplesPerIter:k*numPeriodsPerIter*numSamplePerPeriod)/fs;
    
    if strcmp(robot_movement_type, 'triangular_rastered')
        numPeriod = ceil(FOV_z/time/(driveMag*2)); % number of periods of the triangle wave, this represents the movement in the z direction
        numSamplePerRobotArmPeriod = fs/numPeriod; % sample per period in Fs
        p = 1/numPeriod;

        FFP_x = FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
    elseif strcmp(robot_movement_type, 'linear_rastered')
        x = linspace(-FOV_x/2, FOV_x/2, numLines); % robot arm movement in x direction w.r.t. time
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
        FFP_x = repmat(x(noLine), [1, numPeriodsPerIter*numSamplePerPeriod+1]); % movement of FFP in x direction, for linear x is constant throughout the line
    else
        x = linspace(-FOV_x/2, FOV_x/2, numLines); % robot arm movement in x direction w.r.t. time
        z = repmat(SPIOlocation(1, 2), [1, numPeriodsPerIter*numSamplePerPeriod+1]);
        FFP_x = repmat(x(noLine), [1, numPeriodsPerIter*numSamplePerPeriod+1]); % movement of FFP in x direction, for linear x is constant throughout the line
    end
    drive = driveMag*cos(2*pi*f_drive*t); % drive field movement
    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)

    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFP_speed = sqrt(zDifference.^2 + xDifference.^2);
    FFP_angle = wrapTo360(round(atan2d(xDifference, zDifference), 1));    

    angletUnique = unique(FFP_angle);
    numAngle = length(angletUnique) ;

  
    fprintf('Writing %d of %d, ',nout,signal_size*numIters);
    chunkSize = min(chunk,signal_size*numIters-nout);
    matObj.FFP_x(count,1:signal_size) = FFP_x(1:downsample:end-1);
    matObj.FFP_z(count,1:signal_size) = FFP_z(1:downsample:end-1);
    matObj.z(count, 1:signal_size) = z(1:downsample:end-1);
    matObj.FFP_speed(count,1:signal_size) = FFP_speed(1:downsample:end);
    matObj.FFP_angle(count,1:signal_size) = FFP_angle(1:downsample:end);
    nout = nout + chunkSize;
    
    
    
    endTime(count) = toc;
    fprintf('Iter: %i of %i, Time: %d\n', count, numIters, endTime(count));
    count = count + 1;
end

totalTime = sum(endTime)

FFP_x = transpose(matObj.FFP_x); FFP_x = FFP_x(:);
FFP_z = transpose(matObj.FFP_z); FFP_z = FFP_z(:);
z = transpose(matObj.z); z = z(:);
horizontalSignal_mpi_mat = transpose(matObj.horizontalSignal_mpi_mat); horizontalSignal_mpi_mat = horizontalSignal_mpi_mat(:);

figure; plot(FFP_x, FFP_z); hold on;
plot(FFP_x, z)
xlabel('x-axis (m)'); ylabel('z-axis(cm)');
