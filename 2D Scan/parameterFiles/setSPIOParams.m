function SPIOparams = setSPIOParams(Physicsparams, pos, tau_val)
    SPIOparams = struct;
    
    SPIOparams.diameter = [25]; % (nm)
    SPIOparams.tau = [tau_val, 3e-6]; % (S)
    
    % spio distribution in 2D
    SPIOparams.SPIOdistribution = zeros(1024, 1024, 2); 
    SPIOparams.SPIOdistribution(pos, 512, 1) = 1;
    SPIOparams.SPIOdistribution(200-20:200+20, 256, 2) = 1;
%     SPIOparams.SPIOdistribution(1:40, (225:275), 1) = 1;
%     SPIOparams.SPIOdistribution((1:40)+100, (225:275), 2) = 1;
    % extent of SPIO distribution 
    SPIOparams.image_FOV_x = 0.05; % (m) 
    SPIOparams.image_FOV_z = 0.05; % (m) 
    SPIOparams.dx = SPIOparams.image_FOV_x/size(SPIOparams.SPIOdistribution,1);  % distance between each pixel (m)
    SPIOparams.dz = SPIOparams.image_FOV_z/size(SPIOparams.SPIOdistribution,2);  

    SPIOparams.r_t = cell(1, length(SPIOparams.diameter));
    for k=1:length(SPIOparams.diameter)
        if (SPIOparams.tau(k) == 0)
            SPIOparams.r_t{k} = 1;
        else
            t = (0:1/Physicsparams.fs:SPIOparams.tau(k)*15);
            SPIOparams.r_t{k} = 1/SPIOparams.tau(k)*exp(-t./SPIOparams.tau(k));
        end
    end
end