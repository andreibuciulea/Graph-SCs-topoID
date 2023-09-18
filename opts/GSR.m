function out = GSR(C,reg)

N = size(C,1);

verbose = false;% reg.verbose;
max_iters = 30;%reg.max_iters;
epsilon = 1e-4;%reg.epsilon;
%max_iters = 15;

if verbose
   disp('  -Starting GSR Low Rank optimization...') 
end

for i = 1:max_iters
    cvx_begin quiet
        variable S_hat(N,N) symmetric
        minimize (norm(S_hat(:),1)) %+ lambda/2*square_pos(norm(C*S_hat - S_hat*C, 'fro')))
        subject to
            abs(diag(S_hat)) <= 1e-6;
            S_hat >= 0;
            S_hat*ones(N,1) >= 1;
            norm(C*S_hat - S_hat*C, 'fro') <= epsilon;
    cvx_end
    if strcmp(cvx_status,'Solved') %%|| strcmp(cvx_status,'Failed')  
        if verbose
            disp(['Break after ' num2str(i) ' iterations'])
        end
        break;
    else
        if verbose
            disp('Updating epsilon ...')
        end
        epsilon = epsilon*2;
        %epsilon = epsilon*2;
    end
end
out.H1 = S_hat;
out.H2 = (generate_B(ones(N)).B0'.*kr(S_hat,S_hat))';
out.H2kr = out.H2;
end