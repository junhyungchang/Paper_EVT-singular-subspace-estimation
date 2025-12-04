% Generate data for Gaussian SBM experiment (Supplement)
clear; clc
tic
m = 1.8e3; 
n1 = 3e3;
% m = 1e3;
% n1 = 1e3;
% n2 = 1.2*n1;
n2 = n1;
n=n1;
c = n1/n2;
r = 2;
maxnorms = zeros(m,1);
targetsig = zeros(m,1);
% targetuc = targetsig;
targetc = targetsig;
%%%%%%%%%%%%%Create random test signal matrix%%%%%%%%%%%%%%
U = normrnd(3, 1, [n1,r]);
% U = ones(n1,r);
U = orth(U);
V = normrnd(3,1,[n2,r]);
% V = ones(n2,r);
V = orth(V);
% V =U;
%%%%%%%%%%%%%%%%%%%%%%%SBM signal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = log(n)^(1/2)/sqrt(n);
q = p/log(n);
sig1 = p*(1-p);
sig2 = q*(1-q);
sig = sig1;
s1 = n*0.5; s2 = n-s1;
D = diag([sig1*ones(1,s1), sig2*ones(1,s2)]);
B = [p , q ; q , p];
s = [s1, s2];
M = sbmmean(B, s);
[U,S] = eigs(M,r);
Sdiag = diag(S);
V = U;
sr = S(r,r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n____Progress Bar____\n")
numbar = 20;
for j = 1:m
    E = normrnd(0, 1, [n1,n2]); % independent but not identical case
    E = E*sqrt(D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mhat = M + E;
    [Uhat,Shat,Vhat] = svds(Mhat,r,'largest','MaxIterations', 1000);  
    Shatuc = Shat;
    Shat = diag(Shat);
    Shat = sqrt((Shat.^2 - (1+c)*sig^2*n2 + sqrt(((1+c)*sig^2*n2-Shat.^2).^2 - 4*c*sig^4*n2^2))/2);
    Shat = diag(Shat);
    Hu = Uhat'*U;
    [X,~,Y] = svd(Hu);
    Ru = X*Y';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perturb = Uhat*Ru - U;
    prownorm = tinorm(perturb);
    if mod(j,m/numbar) == 0
        fprintf("1")
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    Ssampc = diag(Shat)'; % corrected
    % Ssamp = diag(Shatuc)'; % uncorrected
    Spluginsig = Sdiag; % Signal singular values
    % Spluginsamp = Ssamp; % uncorrected Sample singular values 
    Spluginsampc = Ssampc; % corrected Sample singular values 
    VDV = V'*D*V;
    lammatsig = S^2\VDV;
    % lammatuc = Shatuc^2\VDV;
    lammatc = Shat^2\VDV;
    %%%%%%%%%%%%%%%%%%%%%%%%%
%     lamssig = sort(2*sig^2./Spluginsig.^2,"descend");
    lamssig = 2*sort(eig(lammatsig),'descend');
    lam1sig = lamssig(1);
    lamremsig = lamssig(2:r);
    lamfactorsig = 1 - lamremsig/lam1sig;
    lamfactorsig = prod(lamfactorsig);
    lamfactorsig = sqrt(lamfactorsig);
    ansig = sqrt(lam1sig)/(2*sqrt(log(2*n)));
    bnsig = sqrt(lam1sig*log(2*n)) - (log(log(2*n))+log(4*pi))*ansig/2;
    targetsig(j) = (prownorm -bnsig)/ansig + log(lamfactorsig);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     lamssampc = sort(2*sig^2./Spluginsampc.^2,"descend");
    lamssampc = 2*sort(eig(lammatc),'descend');
    lam1sampc = lamssampc(1);
    lamremsampc = lamssampc(2:r);
    lamfactorsampc = 1 - lamremsampc/lam1sampc;
    lamfactorsampc = prod(lamfactorsampc);
    lamfactorsampc = sqrt(lamfactorsampc);
    anc = sqrt(lam1sampc)/(2*sqrt(log(2*n)));
    bnc = sqrt(lam1sampc*log(2*n)) - (log(log(2*n))+log(4*pi))*anc/2;
    targetc(j) = (prownorm -bnc)/anc + log(lamfactorsampc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc



save("sbm_norm.mat","targetsig","targetc")


