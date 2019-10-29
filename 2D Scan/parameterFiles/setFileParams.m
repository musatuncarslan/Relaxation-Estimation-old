function [fileObj, signal_size, chunk] = setFileParams(fileName, Physicsparams, MPIparams, SPIOparams, Simparams)


    fileObj = matfile(fileName);
    fileObj.Properties.Writable = true;

    % save simulation parameters
    fileObj.MPIparams = MPIparams;
    fileObj.SPIOparams = SPIOparams;
    fileObj.Physicsparams = Physicsparams;
    fileObj.Simparams = Simparams;

    % create space for the simulation data
    signal_size = Simparams.numSamplesPerIter/Simparams.downsample+2;
    chunk = Simparams.numSamplesPerIter/Simparams.downsample;
    numIters = Simparams.numIters;
    fileObj.horizontalSignal_mpi_mat(numIters,signal_size) = 0;
    fileObj.FFP_x(numIters,signal_size) = 0;
    fileObj.FFP_z(numIters,signal_size) = 0;
    fileObj.FFP_speed(numIters,signal_size) = 0;
    fileObj.FFP_angle(numIters,signal_size) = 0;

end