function [startIter, endIter, samplesPerIter] = linearRasteredIteration(MPIparams, Physicsparams)

    f_drive = MPIparams.f_drive;
    fs = Physicsparams.fs;
    FOV_z = MPIparams.FOV_z;
    
    
    time = MPIparams.time;
    traversedFOVz = MPIparams.traversedFOVz;

    robotSpeed = FOV_z/time; % robot arm movement speed (m/s)
    total_time = (traversedFOVz(2)-traversedFOVz(1))/robotSpeed;

    numPeriods = round(total_time*f_drive); % number of periods on the selected portion of a single line
    div = divisors(numPeriods); % find the divisors of gcd (these numbers divide both total number of samples and number of periods per line)
    numIters = div(div>100 & div < 1000); % use the middle number of divisors as a rule of thumb for number of iterations to solve the whole problem
    if isempty(numIters)
        numIters = div(end);
        msg = sprintf('Low number of iterations, iteration number is chosen as highest possible, %i.\n', numIters);
        fprintf(msg);
    end 


    numIters = numIters(ceil(end/2));
    periodsPerIter = numPeriods/numIters;
    samplesPerIter = periodsPerIter*fs/f_drive;

    startIter = round(traversedFOVz(1)/f_drive/robotSpeed);
    endIter = startIter+numIters-1;
    
end