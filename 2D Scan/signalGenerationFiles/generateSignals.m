function [signal] = generateSignals(colinearIMG, transverseIMG, FFPparams, MPIparams, SPIOparams, Simparams)

    signal = struct;

    FFP_x = FFPparams.FFP_x;
    FFP_z = FFPparams.FFP_z;
    FFP_speed = FFPparams.FFP_speed;
    FFP_angle = FFPparams.FFP_angle;
    FFP_uniqueAngle = FFPparams.FFP_uniqueAngle;
    
    FOV_x = MPIparams.FOV_x;
    FOV_z = MPIparams.FOV_z;
        
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;
    r_t = SPIOparams.r_t;

    % generate colinear and transverse signals using interpolation
    x = (-FOV_x/2:dx:FOV_x/2);
    z = (-FOV_z/2:dz:FOV_z/2);
    numAngle = length(FFP_uniqueAngle);
    numSPIO = length(SPIOparams.diameter);
    colinearSignal = zeros(numSPIO, length(FFP_angle));
    transverseSignal = zeros(numSPIO, length(FFP_angle));

    for l=1:numSPIO
        for k=1:numAngle
            angleIdx = (FFP_angle == FFP_uniqueAngle(k));
            colinearSignal(l, angleIdx) = interp2(x,z,colinearIMG{l, k}, FFP_x(angleIdx), FFP_z(angleIdx), 'spline');
            transverseSignal(l, angleIdx) = interp2(x,z,transverseIMG{l, k}, FFP_x(angleIdx), FFP_z(angleIdx), 'spline');
        end
    end


    % generate horizontal and vertical signals using relaxed colinear and
    % transverse signals.
    horizontalSignal = zeros(size(colinearSignal));
    verticalSignal = zeros(size(colinearSignal));
    for l=1:numSPIO
        for k = 1:numAngle
            angleIdx = (FFP_angle == FFP_uniqueAngle(k));
            horizontalSignal(l, angleIdx) = colinearSignal(l, angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx)) ...
                - transverseSignal(l, angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx));
            verticalSignal(l, angleIdx) = colinearSignal(l, angleIdx).*FFP_speed(angleIdx).*sind(FFP_angle(angleIdx)) ...
                + transverseSignal(l, angleIdx).*FFP_speed(angleIdx).*cosd(FFP_angle(angleIdx));    
        end
    end

    % create relaxed colinear and transverse signals. Crop the signals 
    % appropriately.
    for k=1:numSPIO
        horizontalSignal_temp = conv(horizontalSignal(k,:), r_t{k});
        horizontalSignal(k,:) = horizontalSignal_temp(1:end-length(r_t{k})+1);
        verticalSignal_temp = conv(verticalSignal(k,:), r_t{k});
        verticalSignal(k,:) = verticalSignal_temp(1:end-length(r_t{k})+1);
    end
    
    signal.colinearSignal = colinearSignal(:, 1:Simparams.downsample:end);
    signal.transverseSignal = transverseSignal(:, 1:Simparams.downsample:end);
    signal.horizontalSignal = horizontalSignal(:, 1:Simparams.downsample:end);
    signal.verticalSignal = verticalSignal(:, 1:Simparams.downsample:end);
    
end