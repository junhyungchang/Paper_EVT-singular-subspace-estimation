% Create simulated data for t-distributed noise
clear; clc
tic
m = 1.8e3;
n1 = 2.5e3;
n2 = 1.2*n1;
%%%%%%%
n=n1;
c = n1/n2;
r = 5;
maxnorms = zeros(m,1);
target1 = zeros(m,1);
target2 = zeros(m,1);
target3 = zeros(m,1);
%%%%%%%%%%%%%Create random test signal matrix%%%%%%%%%%%%%%
U = normrnd(3, 1, [n1,r]);
U = orth(U);
V = normrnd(3,1,[n2,r]);
V = orth(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig = 1; % iid case
% Sdiag = linspace(3*sig*sqrt(n1)*log(n1)^(2/3),1.2*sig*sqrt(n1)*log(n1)^(2/3), r);  % weak signal
Sdiag = linspace(3*n1,n1, r);  % strong signal
sr = Sdiag(r);
smin = Sdiag(1);
S = diag(Sdiag);
M = U*S*V';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n____Progress Bar____\n")
numbar = 20;
parfor j = 1:m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nu1 =4;
    E1 = sqrt(sig) * sqrt((nu1-2)/nu1) * trnd(nu1,[n1,n2]);
    nu2 =5;
    E2 = sqrt(sig) * sqrt((nu2-2)/nu2) * trnd(nu2,[n1,n2]);
    nu3 =10;
    E3 = sqrt(sig) * sqrt((nu3-2)/nu3) * trnd(nu3,[n1,n2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mhat1 = M + E1;
    [Uhat1,Shat1,Vhat1] = svds(Mhat1,r,'largest','MaxIterations', 2000);  
    Mhat2 = M + E2;
    [Uhat2,Shat2,Vhat2] = svds(Mhat2,r,'largest','MaxIterations', 2000);  
    Mhat3 = M + E3;
    [Uhat3,Shat3,Vhat3] = svds(Mhat3,r,'largest','MaxIterations', 2000);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hu1 = Uhat1'*U;
    [X,~,Y] = svd(Hu1);
    Ru1 = X*Y';
    Hu2 = Uhat2'*U;
    [X,~,Y] = svd(Hu2);
    Ru2 = X*Y';
    Hu3 = Uhat3'*U;
    [X,~,Y] = svd(Hu3);
    Ru3 = X*Y';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perturb1 = Uhat1*Ru1 - U;
    prownorm1 = tinorm(perturb1);
    perturb2 = Uhat2*Ru2 - U;
    prownorm2 = tinorm(perturb2);
    perturb3 = Uhat3*Ru3 - U;
    prownorm3 = tinorm(perturb3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lamssampc = sort(2*sig^2./Spluginsampc1.^2,"descend");
    lamssampc = sort(2*sig^2./Sdiag.^2,"descend");
%     lamssampc = 2*eig(lammatc);
    lam1sampc = lamssampc(1);
    lamremsampc = lamssampc(2:r);
    lamfactorsampc = 1 - lamremsampc/lam1sampc;
    lamfactorsampc = prod(lamfactorsampc);
    lamfactorsampc = sqrt(lamfactorsampc);
    anc = sqrt(lam1sampc)/(2*sqrt(log(2*n)));
    bnc = sqrt(lam1sampc*log(2*n)) - (log(log(2*n))+log(4*pi))*anc/2;
    target1(j) = (prownorm1 -bnc)/anc + log(lamfactorsampc);
%     %%%%%%%
    % lamssampc = sort(2*sig^2./Spluginsampc2.^2,"descend");
%     lamssampc = 2*eig(lammatc);
    lam1sampc = lamssampc(1);
    lamremsampc = lamssampc(2:r);
    lamfactorsampc = 1 - lamremsampc/lam1sampc;
    lamfactorsampc = prod(lamfactorsampc);
    lamfactorsampc = sqrt(lamfactorsampc);
    anc = sqrt(lam1sampc)/(2*sqrt(log(2*n)));
    bnc = sqrt(lam1sampc*log(2*n)) - (log(log(2*n))+log(4*pi))*anc/2;
    target2(j) = (prownorm2 -bnc)/anc + log(lamfactorsampc);
%     %%%%%%%%%
    % lamssampc = sort(2*sig^2./Spluginsampc3.^2,"descend");
%     lamssampc = 2*eig(lammatc);
    lam1sampc = lamssampc(1);
    lamremsampc = lamssampc(2:r);
    lamfactorsampc = 1 - lamremsampc/lam1sampc;
    lamfactorsampc = prod(lamfactorsampc);
    lamfactorsampc = sqrt(lamfactorsampc);
    anc = sqrt(lam1sampc)/(2*sqrt(log(2*n)));
    bnc = sqrt(lam1sampc*log(2*n)) - (log(log(2*n))+log(4*pi))*anc/2;
    target3(j) = (prownorm3 -bnc)/anc + log(lamfactorsampc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

save("t_dist.mat", "target1", "target2", "target3")