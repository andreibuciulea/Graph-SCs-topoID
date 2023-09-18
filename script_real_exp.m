addpath("opt")
addpath("utils")
addpath("real_dataset")
clear all
%close all
clc
%rng(1)
sig_type = 'vltr';%type of the signal

%models = {'Sam-Cov','VGR-inv','VGR-L21norm','VGR-simp','VGR-sol'};
models = {'GL','GSR','Chepuri','VGR-simp'};

nM = numel(models);

KK_H1 = cell(10,1);
KK_H2 = cell(10,1);
KK_H2_kr = cell(10,1);

filename = 'data/graph_r_20_norm_P.mat';
data = load(filename);
All_X = data.All_X;
All_H1 = data.All_H1;
All_H2 = data.All_H2;

KK = 1;%1:size(All_H1{1},1);
RR = 15:15;
est_err = zeros(numel(RR),numel(KK),3,nM);
fsc_err = zeros(numel(RR),numel(KK),2,nM);
All_est_H1 = {};
All_est_H2 = {};
tic
for rr = 1:numel(RR)
    est_err_r = zeros(numel(KK),3,nM);
    fsc_err_r = zeros(numel(KK),2,nM);
    Xr = squeeze(All_X{rr});
    H1r = squeeze(All_H1{rr});
    H2r = squeeze(All_H2{rr});
    N = size(Xr,2);
    Est_H1r = zeros(numel(KK),N,N,nM);
    Est_H2r = zeros(numel(KK),N,N^2,nM);
    for kk = 1:numel(KK)
        X = squeeze(Xr(kk,:,:));
        H1 = squeeze(H1r(kk,1,:,:));
        H2 = squeeze(H2r(kk,4,:,:));
        [N,M]=size(X);
        disp(['Graph ' num2str(KK(kk)) ' N=' num2str(N)])
        X0 = zeros(size(X));
        for n = 1:N
            X0(n,:) = (X(n,:)-mean(X(n,:)))/std(X(n,:));
        end
        X = X0/max(max(X0));
        %generate volterra graph signals
        sout = sam_cov(X); 
        H1s = sout.H1/max(max(sout.H1));
        H2s = sout.H2/max(max(sout.H2));
        P = X-H1s*X-H2s*kr(X,X);
        norm(X-H1*X-H2*kr(X,X)-P,"fro");
        norm_H1X = norm(H1*X,'fro');
        norm_P = norm(P,'fro');
        norm_H2X = norm(H2*kr(X,X),'fro');
        
        %
        Est_H1 = zeros(N,N,nM);
        Est_H2 = zeros(N,N^2,nM);
        
        alpha = 2;% 
        delta = 3;% 
        beta = 1e-5;% 
        gamma = 17;% 

        rho = 0.4;%GL
        lambda = 1e-8;%VGR_lasso 0
        %lambda = 1e-4;%VGR_lasso 1
        max_iters = 50;
        tol = 1e-6;
        la1 = 1e-3;
        la2 = 1e-3;
        la3 = 1e1;
        alg_prms = struct('alpha',alpha,'beta',beta,'rho',rho,...
                          'lambda',lambda,'max_iters',max_iters,...
                          'delta', delta, 'gamma', gamma,...
                          'la1',la1,'la2',la2,'la3',la3,'tol',tol,'K',4*N);
        e_e = zeros(3,nM);
        e_f = zeros(2,nM);
        for m = 1:nM
            alg_type = models{m};
            estg = estimate_S(X,P,alg_type,alg_prms); %implemet this with the alternating version.
            H1_hat = estg.H1;H1_hat = H1_hat/max(max(H1_hat));
            H2_hat = estg.H2;H2_hat = H2_hat/max(max(H2_hat));
            H2_kr = estg.H2kr;H2_kr = H2_kr/max(max(H2_kr));
            H1 = H1/max(max(H1));
            H2 = H2/max(max(H2));
            nH1 = norm(H1,"fro")^2;
            nH2 = norm(H2,'fro')^2;
            Est_H1(:,:,m) = H1_hat;
            Est_H2(:,:,m) = H2_hat;
            e_e(1,m) = norm(H1-H1_hat,'fro')^2/nH1;
            e_f(1,m) = fscore(H1,H1_hat);
            e_f(2,m) = fscore(H2,H2_hat);
            e_e(2,m) = norm(H2-H2_hat,'fro')^2/nH2;
            e_e(3,m) = norm(H2-H2_kr,'fro')^2/nH2;
        end
        Est_H1r(kk,:,:,:) = Est_H1;
        Est_H2r(kk,:,:,:) = Est_H2;
        est_err_r(kk,:,:) = e_e;
        fsc_err_r(kk,:,:) = e_f;

%         figure()
%         subplot(2,nM+1,1)
%         imagesc(H1)
%         title('True H1')
%         colorbar()
%         subplot(2,nM+1,nM+2)
%         imagesc(H2)
%         colorbar()
%         title('True H2')
%         for k = 1:nM
%             H1_hat = Est_H1(:,:,k);
%             H2_hat = Est_H2(:,:,k);
%             subplot(2,nM+1,k+1)
%             imagesc(estg.H1)
%             title([models{k} ' err H1: ' num2str(e_e(1,k)) ' est fsc: ' num2str(e_f(1,k))])
%             colorbar()
%             subplot(2,nM+1,k+nM+2)
%             imagesc(estg.H2)
%             title([models{k} ' err H2: ' num2str(e_e(2,k)) ' est fsc: ' num2str(e_f(2,k))])
%             colorbar()
%         end
    end
    est_err(rr,:,:,:) = est_err_r;
    fsc_err(rr,:,:,:) = fsc_err_r;
    All_est_H1{rr} = Est_H1r; 
    All_est_H2{rr} = Est_H2r; 

    fn1 = ['VGR_' num2str(RR(rr)) '.mat'];
    save(fn1,"Est_H1r","Est_H2r");
end
toc