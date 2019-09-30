function [FFPparams] = generateFFP(t, MPIparams, Simparams)

    FOV_z = MPIparams.FOV_z;
    FOV_x = MPIparams.FOV_x;
    ffp_type = MPIparams.ffp_type;
    time = MPIparams.time;

    driveMag = MPIparams.driveMag;
    f_drive = MPIparams.f_drive;
    
    numSamplesPerIter = Simparams.numSamplesPerIter;

    if strcmp(ffp_type, 'linear_rastered')
        FFP_x = repmat(MPIparams.ffpPosition(1), [1, numSamplesPerIter+1]); % movement of FFP in x direction, for linear x is constant throughout the line
        z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
    elseif strcmpi(ffp_type, 'fixed')
        FFP_x = repmat(MPIparams.ffpPosition(1), [1, numSamplesPerIter+1]); % movement of FFP in x direction, for linear x is constant throughout the line
        z = repmat(MPIparams.ffpPosition(2), [1, numSamplesPerIter+1]);
    else
        
    end
    
    drive = driveMag*cos(2*pi*f_drive*t); % drive field movement

    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)

    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFP_speed = sqrt(zDifference.^2 + xDifference.^2);
    FFP_angle = wrapTo360(round(atan2d(xDifference, zDifference), 2));    
    
    FFPparams.FFP_x = FFP_x;
    FFPparams.FFP_z = FFP_z;
    FFPparams.FFP_speed = FFP_speed;
    FFPparams.FFP_angle = FFP_angle;
    FFPparams.FFP_uniqueAngle = unique(FFP_angle);
    
end