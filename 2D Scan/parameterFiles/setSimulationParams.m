function [Simparams] = setSimulationParams(MPIparams, Physicsparams)
    Simparams = struct;
    
    ffp_type = MPIparams.ffp_type;
    
    if strcmpi(ffp_type, 'linear_rastered')
        [startIter, endIter, numSamplesPerIter] = linearRasteredIteration(MPIparams, Physicsparams);
        numIters = endIter - startIter;
    elseif strcmpi(ffp_type, 'fixed')
        numIters = MPIparams.cycle;
        total_time = numIters/MPIparams.f_drive;
        numSamplesPerIter = Physicsparams.fs*total_time/numIters;
        startIter = 1;
        endIter = startIter + numIters - 1;
    elseif strcmpi(ffp_type, 'triangular')
    end
    
    Simparams.startIter = startIter;
    Simparams.endIter = endIter-1;
    Simparams.numSamplesPerIter = numSamplesPerIter;
    Simparams.numIters = numIters;

    Simparams.downsample = Physicsparams.fs/MPIparams.fs; % downsample ratio

end