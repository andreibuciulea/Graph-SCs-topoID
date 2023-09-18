function out = gen_graph_volt_sig(prms)
    %params related to the graph generation
    N = prms.N;%number of nodes
    p = prms.p;%link probability

    %params related to the signal generation
    M = prms.M;%number of samples
    T = 20;%max number of iterations for data generation
    K = 200;
    tol_H1 = 0.3;
    tol_x = 1e-8;
    
    %find a valid Adjacency matrix fot data generation
    valid = false;
    i = 1;
    while ~valid
        H1 = generate_connected_ER(N,p);
        valid = min(abs(eig(H1)-ones(N,1))) > tol_H1;
        i = i+1;
        if i > K
            error('No valid adjacency, max tries reached')
        end
    end
    
    %Generate second order interacion matrix 
    Bout = generate_B(H1);
    %H2 = Bout.B2;
    H2 = Bout.B1;
    %Additional variables (initialization)
    m = 10;
    R = ceil(M/m);
    In = eye(N);Im = eye(m);Imn = eye(m*N);On = ones(1,N);
    B = kron(Im,In-H1);
    B_inv = kron(Im,inv(In-H1));
    %x1_x2_norm = zeros(T,R,2);
    X1M = zeros(N,R*m);
    PM = zeros(N,R*m);
    for r = 1:R
        %generate exogenous data and initialization
        P = 0.01*randn(N,m);
        X1 = randn(N,m);
        
        %Generate data
        for t = 1:T
            C1 = kron(Im,H2)*kr(Imn,X1*kron(Im,On));
            C2 = kron(Im,H2)*kron(kr(Im,X1),In);
            %compute x1 using C1 and C2
            F = B_inv*C2/(B-C1)+ B_inv;
            %x1 = R*vec(P);
            %reshape x1 into a matrix
            X2 = reshape(F*vec(P),[N,m]);
            a = norm(X1-X2,"fro");
            %x1_x2_norm(t,r,1) = a;
            %x1_x2_norm(t,r,2) = norm(X1-H1*X1-H2*kr(X1,X1)-P,"fro");
            X1 = X2;
    
            if a <= tol_x
                %disp(['Tolerance achieved at ' num2str(t) ' iteration'])
                break;
            end
    
            if t == T && a > tol_x
                disp('Invalid adjacency matrix for data generation!')
                out = gen_graph_volt_sig(prms);
                break;
            end
            X1M(:,(r-1)*m+1:r*m) = X1;
            PM(:,(r-1)*m+1:r*m) = P;
        end
    end
    %X1M = randn(N,M);
    %PM = X1M-H1*X1M-H2*kr(X1M,X1M);
    out.P = PM;
    out.X = X1M;
    out.H1 = H1;
    out.H2 = H2;
    %out.x1_x2_norm = x1_x2_norm;
end