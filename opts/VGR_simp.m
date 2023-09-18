function out = VGR_simp(X,P,prms)
    N = size(X,1);
    alpha = 5;%prms.alpha;
    beta = 1e-5;%prms.beta;
    max_iters = prms.max_iters;
    B0 = generate_B(ones(N)).B0;
    H2 = B0;
    v = zeros(max_iters,3);
    In = eye(N);
    la = 1e-1;
    tol = 1e-10;
    for t = 1:max_iters
   
    % Estimate H1
        C1 = X-H2*kr(X,X)-P; 
        h1 = kron(inv(X*X'+In),In)*(vec(C1*X')-alpha/2*ones(N^2,1));%for real data estimations
        %h1 = kron(inv(X*X'),In)*(vec(C1*X')-alpha/2*ones(N^2,1));
        H1 = reshape(h1,[N,N]);
    
    % Proyect H1 into valid adjacency matrices
        H1 = H1 - diag(diag(H1));
        H1 = (H1+H1')/2;
        H1(H1<0) = 0;
        %H1 = H1/max(max(H1));
    
    %Estimate H2
        C2 = X-H1*X-P;
        Y = kr(X,X);
        h2 = kron(inv(Y*Y'+la*eye(N^2)),In)*(vec(C2*Y')-beta/2*ones(N^3,1));
        H2 = reshape(h2,[N,N^2]);
        
    %Proyect H2 into a valid second-order relationship (triangles) matrix 
        R = 1:N;
        for k = 1:N
            for i = R(R~=k)
                for j = R(R~=k & R~=i)
                    %this can be solved more efficiently, by assuming symmetry. 
                    H2(k,N*(i-1)+j) = H1(k,i)*H1(k,j)*H1(j,i);
                    %H2(k,N*(i-1)+j) = H1(k,i)*H1(k,j);
                end
            end
        end
        H2 = H2.*B0;
        H2(H2<0) = 0;
        %H2 = H2/max(max(H2));

        v(t,1) = norm(X - H1*X-H2*kr(X,X)-P,'fro');
        if t > 1
            if abs(v(t-1,1)-v(t,1)) < tol
                break;
            end
        end
        v(t,2) = norm(H1,1);
        v(t,3) = norm(H2,1);
    end
    %figure(2);semilogy(v);legend('1','2','3');grid on;
    out.H1 = H1;
    out.H2 = H2; 
    out.H2kr = (B0'.*kr(H1,H1))';
end