function [colinearPSF,transversePSF, X, Z] = generatePSF(MPIparams, SPIOparams, FFPparams)

    G = [MPIparams.Gxx, MPIparams.Gyy, MPIparams.Gzz];
    diameter = SPIOparams.diameter;
    FOV_x = MPIparams.FOV_x; 
    FOV_z = MPIparams.FOV_z;
    dx = SPIOparams.dx; 
    dz = SPIOparams.dz;

    FFP_uniqueAngle = FFPparams.FFP_angle;
    numAngle = length(FFPparams.FFP_uniqueAngle);

    x = (-FOV_x/2:dx:FOV_x/2);
    z = (-FOV_z/2:dz:FOV_z/2);
    [X,Z] = meshgrid(x,z);
    Y = 0;

    mo = 4*pi*1e-7;
    Gxx = G(1)/mo;
    Gyy = G(2)/mo;
    Gzz = G(3)/mo;
    
    colinearPSF = cell([length(diameter), numAngle]);
    transversePSF = cell([length(diameter), numAngle]);
    
    for m=1:length(diameter)
        d = diameter(m)*1e-9;
        kb = 1.3806488e-23;
        T = 300;
        k = (0.1*pi*d^3)/(kb*T);

        for l=1:numAngle
            e = transpose([X(:), zeros(size(X,1)*size(X,2), 1), Z(:)]);
            R = roty(FFP_uniqueAngle(l));

            e = R*e;
            rotX = reshape(e(1, :), size(X,1), size(X,2));
            rotZ = reshape(e(3, :), size(Z,1), size(Z,2));

            Hxyz = sqrt((Gxx*rotX).^2 + (Gyy*Y).^2 + (Gzz*rotZ).^2);
            Hxyzk = Hxyz*k;
            idx = (Hxyz == 0);


            Lx = zeros(size(Hxyz));
            Lx(~idx) = coth(Hxyzk(~idx)) - 1./(Hxyzk(~idx)) ;
            Lx(idx) = 0;

            Lx_der = zeros(size(Hxyz));
            Lx_der(~idx) = 1./Hxyzk(~idx).^2 - 1./sinh(Hxyzk(~idx)).^2;
            Lx_der(idx) = 1/3;


            Tenv = Lx_der; % tangential component
            Nenv = Lx./Hxyzk; Nenv(isnan(Nenv)) = 1/3; % normal component
            arg_colli = Gzz^2.*rotZ.^2./Hxyz.^2; arg_colli(isnan(arg_colli)) = 0;

            colinearPSF{m, l} = Tenv.*arg_colli + Nenv.*(1-arg_colli);

            arg_trans = Gxx*Gzz*rotX.*rotZ./Hxyz.^2; arg_trans(isnan(arg_trans)) = 0;
            transversePSF{m, l} = (Tenv - Nenv).*arg_trans;
        end
    end
end

