function out = estimate_S(X,P,model,prms)
    switch model
        case 'VGR-inv'
            out = VGR_inv(X,P);
        case 'Sam-Cov'
            out = sam_cov(X);
        case 'GL'
            C = (X-P)*(X-P)';
            %C = X*X';
            %C = C/max(max(C));
            %C = inv(sqrt(diag(diag(C))))*C*inv(sqrt(diag(diag(C))));
            out = graphicalLasso(C,prms);
        case 'GSR'
            C = (X-P)*(X-P)';
            %C = inv(sqrt(diag(diag(C))))*C*inv(sqrt(diag(diag(C))));
            C = C/max(max(C));
            out = GSR(C,prms);
        case 'VGR-lasso-0'
            out = VGR_lasso(X,P,prms,0);
        case 'VGR-lasso-1'
            out = VGR_lasso(X,P,prms,1);
        case 'VGR-alt'
            out = VGR_alt(X,P,prms);
        case 'VGR-alt-rw'
            out = VGR_alt_rw(X,P,prms);
        case 'VGR-simp'
            out = VGR_simp(X,P,prms);
        case 'VGR-L21norm'
            out = VGR_L21norm(X,P,prms);
        case 'VGR-sol'
            out = VGR_sol(X,P,prms);
        case 'VGR-it'
            out = VGR_it(X,P,prms);
        case 'VGR-group'
            out = VGR_group(X,P,prms);
        case 'VGR-group-v1'
            out = VGR_group_v1(X,P,prms);
        case 'VGR-2group'
            out = VGR_2group(X,P,prms);
        case 'VGR-fast'
            out = VGR_fast(X,P,prms);
        case 'VGR-triang'
            out = VGR_triang(X,P,prms);
        case 'Chepuri'
            out = chepuri_alg(X,prms);
        otherwise
            disp('Unknown Method')
    end
end
