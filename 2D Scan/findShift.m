function del_t = findShift(R, f_d, MPIparams)
    arg = 2*pi*f_d;

    % taylor series expension of sine
    N = 11;
    p = zeros(1, N+1); 
    for k=0:N
        p(k+1) = ((-1)^k)/factorial(2*k+1)*(arg^(2*k+1)); 
    end
    p = p*MPIparams.driveMag;
    p = flip(p);
    p = upsample(p, 2); % turn into polynomial
    p(end-1) = p(end-1) - R; % add the extra R*del_t 
    rp = roots(p); % find roots
    real_idx = imag(rp) == 0; % extract real roots
    real_roots = rp(real_idx);
    del_t = real_roots(end); % select the last one (last one is always the root we are looking for)

end