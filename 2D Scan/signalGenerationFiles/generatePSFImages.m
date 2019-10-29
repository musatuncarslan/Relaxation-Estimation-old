function [colinearIMG, transverseIMG] = generatePSFImages(SPIOparams, colinearPSF, transversePSF)

    numAngle = size(colinearPSF, 2);
    numParticle = size(colinearPSF, 1);

    colinearIMG = cell(numParticle, numAngle);
    transverseIMG = cell(numParticle, numAngle);
    for m=1:numParticle
        tempDistribution = SPIOparams.SPIOdistribution(:,:,m);

        % pad the distribution with 0s so it has same size with the PSFs
        img_size = size(tempDistribution); psf_size = size(colinearPSF{1});
        if img_size(1) < psf_size(1)
           tempDistribution = padarray(tempDistribution,[floor((psf_size(1)-img_size(1))/2), 0],0,'both');
           img_size = size(tempDistribution);
           tempDistribution = padarray(tempDistribution,[psf_size(1)-img_size(1), 0],0,'post');
        end
        if img_size(2) < psf_size(2)
           tempDistribution = padarray(tempDistribution,[0, floor((psf_size(2)-img_size(2))/2)],0,'both');
           img_size = size(tempDistribution);
           tempDistribution = padarray(tempDistribution,[0, psf_size(2)-img_size(2)],0,'post');
        end
        % multiplication in frequency domain is faster than convolution in
        % time domain
        imgF = fft2(tempDistribution, size(colinearPSF{1},1), size(colinearPSF{1},2));
        for k=1:numAngle
            colinearIMG{m,k} = fftshift(ifft2(imgF.*fft2(colinearPSF{m,k})));    
            transverseIMG{m,k} = fftshift(ifft2(imgF.*fft2(transversePSF{m,k})));   
        end
    end
end