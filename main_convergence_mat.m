% Creates simulated data for the three histogram-qqplot pairs in Figure 1.
clear; clc
tic
m = 1.8e3; 
n1 = 3e3;
% m = 5e2; % smaller test
% n1 = 7e2;
n2 = 1.2*n1;
n=n1;
c = min([n1/n2, n2/n1]);
r = 5;
p = 1;
maxnorms = zeros(m,1);
targetsig = zeros(m,1);
targetuc = targetsig;
targetc = targetsig;
%%%%%%%%%%%%%Create signal matrix%%%%%%%%%%%%%%
U = normrnd(3, 1, [n1,r]);
U = orth(U);
V = normrnd(3,1,[n2,r]);
V = orth(V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig = 1; % iid case
% Sdiag = linspace(6*sqrt(n1),4*sqrt(n1), r);  % tall
% Sdiag = linspace(7*sqrt(n2),5*sqrt(n2), r);  % wide
% Sdiag =
% % % % % weak signal 1 (old normalizing sequences)
% linspace(2*sig*sqrt(n1)*log(n1)^(1.1),0.5*sig*sqrt(n1)*log(n1)^(1.1), r);
% % % % % weak signal 2 (new)
% Sdiag = linspace(1*sig*sqrt(r*n1)*log(n1)^(1.01),0.25*sig*sqrt(r*n1)*log(n1)^(1.01), r);  
% % % % % p-repeated weak sr
% Sdiag = [3*sig*sqrt(r*n1)*log(n1)^(1.01)*ones(1, r-p),0.25*sig*sqrt(r*n1)*log(n1)^(1.01)*ones(1,p)]; 
% % % % % strong signal
Sdiag = linspace(3*n1,n1, r); 
% % % % % p-repeated strong sr
% Sdiag = [linspace(2*n1, 3*n1, r-p), n1*ones(1,p)];  
sr = Sdiag(r);
smin = Sdiag(1);
S = diag(Sdiag);
M = U*S*V';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n____Progress Bar____\n")
numbar = 20;
for j = 1:m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E = normrnd(0, sig, [n1,n2]); % iid Gaussian case
    Mhat = M + E;
    [Uhat,Shat,Vhat] = svds(Mhat,r,'largest','MaxIterations', 1000);  
    Shatuc = Shat;
    Shat = diag(Shat);
    N = max([n1,n2]);
    Shat = sqrt((Shat.^2 - (1+c)*sig^2*N + sqrt(((1+c)*sig^2*N-Shat.^2).^2 - 4*c*sig^4*N^2))/2);
    Shat = diag(Shat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hu = Uhat'*U;
    [X,~,Y] = svd(Hu);
    Ru = X*Y';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    perturb = Uhat*Ru - U;
    prownorm = tinorm(perturb);
    approx = E*V/S;
    if mod(j,m/numbar) == 0
        fprintf("1")
    end
    %%%%%%%%%%%iid entries case%%%%%%%%%%%%
    Ssampc = diag(Shat)'; % corrected
    Ssamp = diag(Shatuc)'; % uncorrected
    Spluginsig = Sdiag; % Signal singular values
    Spluginsamp = Ssamp; % Sample singular values
    Spluginsampc = Ssampc;
    %%%%%%%%%%%%%%Oracle%%%%%%%%%%%%%%%%%%
    lamssig = sort(2*sig^2./Spluginsig.^2,"descend"); % iid entries
    lam1sig = lamssig(1);
    lamremsig = lamssig(p+1:r);
    lamfactorsig = 1 - lamremsig/lam1sig;
    lamfactorsig = prod(lamfactorsig);
    lamfactorsig = sqrt(lamfactorsig);
    ansig = sqrt(lam1sig)/(2*sqrt(log(n)));
    bnsig = sqrt(lam1sig * log(n)) ...
        + (p-2)*sqrt(lam1sig) *log(log(n))/(4*sqrt(log(n))) ...
        - sqrt(lam1sig)*gammaln(p/2)/(2*sqrt(log(n)));
    targetsig(j) = (prownorm -bnsig)/ansig + log(lamfactorsig);
    %%%%%%%%%%%%Uncorrected%%%%%%%%%%%%%%%%%%%%%
    lamssamp = sort(2*sig^2./Spluginsamp.^2,"descend"); % iid entries
    lam1samp = lamssamp(1);
    lamremsamp = lamssamp(p+1:r);
    lamfactorsamp = 1 - lamremsamp/lam1samp;
    lamfactorsamp = prod(lamfactorsamp);
    lamfactorsamp = sqrt(lamfactorsamp);
    anuc = sqrt(lam1samp)/(2*sqrt(log(n)));
    bnuc = sqrt(lam1samp * log(n)) ...
        + (p-2)*sqrt(lam1samp) *log(log(n))/(4*sqrt(log(n))) ...
        - sqrt(lam1samp)*gammaln(p/2)/(2*sqrt(log(n)));
    targetuc(j) = (prownorm -bnuc)/anuc + log(lamfactorsamp);
    %%%%%%%%%%%Corrected%%%%%%%%%%%%%%%%%%%%
    lamssampc = sort(2*sig^2./Spluginsampc.^2,"descend"); % iid entries
    lam1sampc = lamssampc(1);
    lamremsampc = lamssampc(p+1:r);
    lamfactorsampc = 1 - lamremsampc/lam1sampc;
    lamfactorsampc = prod(lamfactorsampc);
    lamfactorsampc = sqrt(lamfactorsampc);
    anc = sqrt(lam1sampc)/(2*sqrt(log(n)));
    bnc = sqrt(lam1sampc * log(n)) ...
        + (p-2)*sqrt(lam1sampc) *log(log(n))/(4*sqrt(log(n))) ...
        - sqrt(lam1sampc)*gammaln(p/2)/(2*sqrt(log(n)));
    targetc(j) = (prownorm -bnc)/anc + log(lamfactorsampc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc

save("main_convergence.mat","targetsig","targetc","targetuc")
