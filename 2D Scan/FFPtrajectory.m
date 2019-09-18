function [FFP_x, FFP_z, FFPspeed, FFPangle, x, z] = FFPtrajectory(FOV_x, FOV_z, fs, f_drive, time, driveMag, period)
    

    numPeriod = ceil(FOV_z/time/(driveMag*2)); % number of periods of the triangle wave, this represents the movement in the z direction
    numSamplePerRobotArmPeriod = fs/numPeriod; % sample per period in Fs
    periodIdx = period*numSamplePerRobotArmPeriod; % period index of concern
    p = 1/numPeriod;
    
    t = ((periodIdx(1):periodIdx(end))/fs);

    x = FOV_x*(2*abs(2*(t/p - floor(t/p + 0.5)))-1)/2; % robot arm movement in x direction w.r.t. time
    z = FOV_z/time*t - FOV_z/2; % robot arm movement in z direction w.r.t. time
    drive = driveMag*cos(2*pi*f_drive*t); % drive field movement

    FFP_x = x; % movement of FFP in x direction
    FFP_z = z+drive; % movement of FFP in z direction (robot arm movement in z direction + drive field movement, since drive field is in z direction as well)

    zDifference = diff(FFP_z, 1, 2);
    xDifference = diff(FFP_x, 1, 2);

    FFPspeed = sqrt(zDifference.^2 + xDifference.^2);
    FFPangle = wrapTo360(round(atan2d(xDifference, zDifference), 1));

   
    

end