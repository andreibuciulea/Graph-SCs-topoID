function out = generate_B(A)
    N = size(A,1);
    B0 = ones(N,N^2);
    B0(:,1:N+1:N^2) = 0;
    for k = 1:N
        B0(k,N*(k-1)+1:N*k) = 0;
        B0(k,k:N:N^2) = 0;
    end
%     noB = ~B;
%     figure()
%     subplot(121)
%     imagesc(B)
%     title('Matrix B')
%     subplot(122)
%     imagesc(noB)
%     title("Matrix no B")
    B1 = B0;
    B2 = B0;
    
    %To have a triangle, all interactions between the 3 nodes are needed
    R = 1:N;
    for k = 1:N
        for i = R(R~=k)
            for j = R(R~=k & R~=i)
                %this can be solved more efficiently, by assuming symmetry. 
                B1(k,N*(i-1)+j) = A(k,i)*A(k,j)*A(j,i);
                B2(k,N*(i-1)+j) = A(k,i)*A(k,j);
            end
        end
    end
    out.B0 = B0;
    out.B1 = B1;
    out.B2 = B2;

%     figure(1)
%     subplot(121)
%     imagesc(B)
%     title('B assuming triangles')
%     subplot(122)
%     imagesc(B2)
%     title('B assuming more interacions')
%     grid on
end

