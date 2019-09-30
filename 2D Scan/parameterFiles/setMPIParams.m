function MPIparams = setMPIParams(Physicsparams)
    MPIparams = struct;


    MPIparams.ffp_type = 'fixed';
    % gradients (T/m) (current scanner)
    MPIparams.Gxx = 4.8;
    MPIparams.Gyy = 2.4;
    MPIparams.Gzz = 2.4;
    
    MPIparams.Bp = 15e-3; % Drive field (T)
    MPIparams.f_drive = 10e3; % drive field frequency

    Hp=MPIparams.Bp/Physicsparams.mu0; % magnetization moment
    G=MPIparams.Gzz/Physicsparams.mu0; % gradient
    MPIparams.driveMag=Hp/G; % extent of the drive field

    MPIparams.fs = 2e6; % sample frequency of the MPI system (Hz)
    MPIparams.FOV_z = 0.06; % FOV in z-axis (meters) (bore axis)
    MPIparams.FOV_x = 0.05; % FOV in x-axis (meters)
    
    % for linear rastered
    MPIparams.time = 0.006; % time (seconds)
    MPIparams.traversedFOVz = [-0.03 0.03]; % traversed fov in the simulation in z-axis (m)
    % for fixed
    MPIparams.cycle = 5; % number of cycles on the fixed position
    MPIparams.ffpPosition = [0, 0]; % ffp position in x and z coordinates
   
end