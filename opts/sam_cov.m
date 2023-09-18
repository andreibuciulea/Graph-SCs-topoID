function out = sam_cov(X)
    N = size(X,1);
    S = X*X';
    S = S-diag(diag(S));
    S(S<0) = 0;
    out.H1 = S;
    out.H2 = (generate_B(ones(N)).B0'.*kr(S,S))';
    out.H2kr = out.H2;
end