function [Simparams] = setSimulationParams(MPIparams, Physicsparams)
    Simparams = struct;
    
    ffp_type = MPIparams.ffp_type;
    
    if strcmpi(ffp_type, 'linear_rastered')
        [startIter, endIter, samplesPerIter] = linearRasteredIteration(MPIparams, Physicsparams);
        numIters = endIter - startIter+1;
        Simparams.startIter = startIter;
        Simparams.endIter = endIter;
    elseif strcmpi(ffp_type, 'fixed')
        numIters = MPIparams.cycle;
        total_time = numIters/MPIparams.f_drive;
        samplesPerIter = Physicsparams.fs*total_time/numIters;
        Simparams.startIter = 1;
        Simparams.endIter = Simparams.startIter + numIters - 1;
    elseif strcmpi(ffp_type, 'triangular')
    end
    

    Simparams.numSamplesPerIter = samplesPerIter;
    Simparams.numIters = numIters;

    Simparams.downsample = Physicsparams.fs/MPIparams.fs; % downsample ratio

end