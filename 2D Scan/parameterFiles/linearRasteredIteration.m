function [startIter, endIter, numSamplesPerIter] = linearRasteredIteration(MPIparams, Physicsparams)

    f_drive = MPIparams.f_drive;
    fs = Physicsparams.fs;
    FOV_z = MPIparams.FOV_z;
    
    
    time = MPIparams.time;
    traversedFOVz = MPIparams.traversedFOVz;


    f_drive_temp = f_drive;
    divider = 1;
    count = -1;
    while (mod(f_drive_temp, divider) == 0)
        count = count + 1;
        divider = divider*10;   
    end

    robotSpeed = FOV_z/time; % robot arm movement speed (m/s)
    bFOVz = -FOV_z/2; % beginning point of FOV in z-axis (m)
    t1 = round((traversedFOVz(1)-bFOVz)/robotSpeed, count);
    t2 = round((traversedFOVz(2)-bFOVz)/robotSpeed, count);
    total_time = round(t2-t1, count);

    numPeriod = total_time*f_drive; % number of periods on the selected portion of a single line
    a = gcd(fs*total_time, numPeriod); % find the greatest common divisor of total number of samples and total number of periods per line
    div = divisors(a); % find the divisors of gcd (these numbers divide both total number of samples and number of periods per line)
    numIters = div(div>100 & div < 1000); % use the middle number of divisors as a rule of thumb for number of iterations to solve the whole problem
    if isempty(numIters)
        numIters = div(end);
        msg = sprintf('Low number of iterations, iteration number is chosen as highest possible, %i.\n', numIters);
        fprintf(msg);
    end 


    numIters = numIters(ceil(end/2));
    numSamplesPerIter = fs*total_time/numIters;

    startIter = (t1*fs)/numSamplesPerIter;
    endIter = startIter+numIters;
    
end