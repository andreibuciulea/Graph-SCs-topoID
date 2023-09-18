addpath("opt")
addpath("utils")
clear all
%close all
clc
rng(7)
nG = 1;
Samp = [50,100,200,300,400,500];
nS = numel(Samp);
N = 20; %number of nodes
sig_type = 'vltr';%type of the signal
%models = {'Sam-Cov','VGR-inv','GL','GSR','VGR-2group'};%algoritm type for graph estimation
%models = {'VGR-simp','VGR-group','VGR-group-v1','VGR-triang'};
models = {'GL','GSR','Chepuri','VGR-2group'};
nM = numel(models);

p = 0.3;%link probability for ER graphs
%params for generate graph and volterra graph signals
prms=struct('N',N,'M',Samp(end),'p',p);%parameters
%amount of noise
sigma = 0;

K = 55; %Number of edges for chepuris algorithm
%
%all_H1 = zeros(N,N,nM);
%all_H2 = zeros(N,N^2,nM);
est_err = zeros(nG,nS,nM,3);
tic
%estimate the graph topology
alpha = 1e-8;%H1 sparsity
beta = 1e-6;%row sparsity
delta = 0;%H2 sparsity
gamma = 1e-6;%column sparsity
rho = 1e-5;%GL
%lambda = 1e-8;%VGR_lasso 0
lambda = 1e-4;%VGR_lasso 1
max_iters = 10;
tol = 1e-6;
la1 = 1e-3;
la2 = 1e-3;
la3 = 1e-3;
alg_prms = struct('alpha',alpha,'beta',beta,'rho',rho,...
                  'lambda',lambda,'max_iters',max_iters,...
                  'delta',delta,'gamma',gamma,...
                  'la1',la1,'la2',la2,'la3',la3,'tol',tol,'K',K);
%load('input_data_exp1.mat','all_GS');
for g = 1:nG
    est_err_g = zeros(nS,nM,3);
    %generate signals
    GS = gen_graph_volt_sig(prms);
    %GS = all_GS{g};
    X0 = GS.X;P0 = GS.P;H1 = GS.H1;H2 = GS.H2;
    %add noise
    AN = add_noise(X0,sigma);
    X0 = AN.Xn;
    g
    for s = 1:nS
        X = X0(:,1:Samp(s));
        P = P0(:,1:Samp(s));
        for m = 1:nM
            alg_type = models{m};
            estg = estimate_S(X,P,alg_type,alg_prms); %implemet this with the alternating version.
            H1_hat = estg.H1;H1_hat = H1_hat/max(max(H1_hat));
            H2_hat = estg.H2;H2_hat = H2_hat/max(max(H2_hat));
            H2_kr = estg.H2kr;H2_kr = H2_kr/max(max(H2_kr));
            nH1 = norm(H1,"fro")^2;
            nH2 = norm(H2,'fro')^2;
            est_err_g(s,m,1) = norm(H1-H1_hat,'fro')^2/nH1;
            est_err_g(s,m,2) = norm(H2-H2_hat,'fro')^2/nH2;
            est_err_g(s,m,3) = norm(H2-H2_kr,'fro')^2/nH2;
%             figure();
%             sgtitle(alg_type)
%             subplot(221)
%             imagesc(estg.H1)
%             colorbar()
%             subplot(222)
%             imagesc(estg.H2)
%             colorbar()
%             subplot(223)
%             imagesc(H1)
%             colorbar()
%             subplot(224)
%             imagesc(H2)
%             colorbar()
        end

    end
    est_err(g,:,:,:) = est_err_g;
end
toc