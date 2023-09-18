function out = VGR_2group(X,P,prms)
    N = size(X,1);
    alpha = prms.alpha;
    delta = prms.delta;
    beta = prms.beta;
    gamma = prms.gamma;
    M = [X; kr(X,X)]';
    Bout = generate_B(ones(N));
    B0 = ~Bout.B0;
    cvx_begin quiet
        variable theta(N+N^2,N)
        F = 0;
        for n = 1:N
            xn = (X(n,:)-P(n,:))';
            F = F + square_pos(norm(xn - M*theta(:,n),2)) + alpha*norm(theta(:,n),1)...
                + beta*sum(norms(reshape(theta(:,n),[N,N+1])',2));
        end
        F = F + gamma*sum(norms([vec(theta(1:N,1:N)) theta(N+1:end,:)]',2))...
            + 0*norm(theta(1:N,:),1) + delta*norm(theta(N+1:end),1); 
        minimize(F)
        subject to 
            theta >= 0;
            diag(theta(1:N,:)) <= 1e-6;%zeros in the ajacency matrix
            %theta(1:N,:)==(theta(1:N,:)+theta(1:N,:)')/2;
            for n = 1:N
                 theta(N*(n-1)+1:n*N,:) == theta(N*(n-1)+1:n*N,:)'; 
            end
            
            theta(N+1:end,:)'.*B0 <= 1e-6;
    cvx_end
    H1 = theta(1:N,:)';
    out.H1 = H1;
    out.H2 = theta(N+1:end,:)'; 
    out.H2kr = (Bout.B0'.*kr(H1,H1))';
end

