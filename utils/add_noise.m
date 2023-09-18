function out = add_noise(X,sigma)
    [N,M] = size(X);
    p_x = norm(X, 'fro')^2/M;
    sigma = sqrt(sigma*p_x/N);
    noise = randn(N, M)*sigma;
    Xn = X + noise;
    p_n = norm(noise, 'fro')^2/M;
    snr = p_x/p_n;

    out.noise = noise;
    out.snr = snr;
    out.Xn = Xn;
    out.p_x= p_x;
    out.p_n = p_n;
end