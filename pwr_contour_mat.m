% % Create data matrix for power contour (Figure 3)
% % $n$ ranges from 20 to 400 and $d_{n}$ varies.
clear; clc
tic
m = 4e2; % choose value smaller than n1
numbar = 20;
reps = 5;
ns = linspace(20,400,20);
numds = ceil(100./ns.^(1/3));
r = 10;
alpha= 0.05;
qgum = gumbel_quantile(1-alpha);
powermat = zeros(sum(numds), 3);
index = 0;
for n = ns
    n1 = n;
    n2 = 1.2*n1;
    c = n1/n2;
    index = index +1;
    %%%%%%%%%%%%%Create random test signal matrix%%%%%%%%%%%%%%
    U0 = orth(normrnd(2,5,[n1,r]));
    V = orth(normrnd(0,5,[n2,r]));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig = 1; % iid case
%     Sdiag = linspace(3*sig*sqrt(n1)*log(n1)^(2/3),1.2*sig*sqrt(n1)*log(n1)^(2/3), r);  % weak signal
    Sdiag = linspace(3*n1^(1),n1^(1), r);  % strong signal
    sr = Sdiag(r);
    smin = Sdiag(1);
    S = diag(Sdiag);
    ts = linspace(0,1,numds(index));
    for l = 1:numds(index)
        [U,d] = generate_U1(U0, ts(l));
        M = U*S*V';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf("\n n = %.d (%.d / %.d), numd = %.d / %.d", n1, index, length(numds), l, numds(index))
        fprintf("\n____Progress Bar____\n")
        powerour = zeros(reps,1);
        for k = 1:reps
            tstats = zeros(m,1);
            parfor j = 1:m
                E = normrnd(0, sig, [n1,n2]); % iid case
                Mhat = M + E;
                [Uhat,Shat,Vhat] = svds(Mhat,r,'largest','MaxIterations', 1000); 
                Shat = diag(Shat);
                Shat = sqrt((Shat.^2 - (1+c)*sig^2*n2 + sqrt(((1+c)*sig^2*n2-Shat.^2).^2 - 4*c*sig^4*n2^2))/2);
                Shat = diag(Shat);
                %%%%%%%%%%%%Our test%%%%%%%%%%%%%%%%%%%
                Hu = Uhat'*U0;
                [X,~,Y] = svd(Hu);
                Ru = X*Y';
                target = tinorm(Uhat*Ru - U0);
                Splugin = diag(Shat);
                lams = sort(2*sig^2./Splugin.^2,"descend");
                lam1 = lams(1);
                lamrem = lams(2:r);
                lamfactor = 1 - lamrem/lam1;
                lamfactor = prod(lamfactor);
                lamfactor = sqrt(lamfactor);
                an = sqrt(lam1)/(2*sqrt(log(2*n)));
                bn = sqrt(lam1*log(2*n)) - (log(log(2*n))+log(4*pi))*an/2 ;
                tstats(j) = (target -bn)/an + log(lamfactor);
            end
            powerour(k) = length(tstats(tstats>=qgum))/m;
            fprintf("1111")
        end 
        powermat(sum(numds(1:index-1))+l,:) = [n1, d, mean(powerour)];
    end
end
toc

buffer = [ns', 1.5 * ones(length(ns),1), nan(length(ns),1)];
powermat = [powermat; buffer];

save("pwr_contour.mat","powermat")

