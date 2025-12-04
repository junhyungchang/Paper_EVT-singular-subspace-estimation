% Create simulated data for power comparison contour plot (Figure 2)
clear; clc
tic
m = 6e2; % number of MC samples
reps = 5;
n1 = 4e2;
n2 = 1*n1;
c = n1/n2;
n= n1;
nummus = 11;
mus = ceil(linspace(380,n1,nummus));
mus(nummus) = n1-1;
ss = n.^(linspace(4,10,9)*0.1);
r = 1;
alpha= 0.05;
qnorm = norminv(1-alpha);
qgum = gumbel_quantile(1-alpha);
powerfrobmat = zeros(length(mus)*length(ss),3);
powerourmat = powerfrobmat;
indexmu = 0;
for mu = mus
    indexmu = indexmu + 1;
    %%%%%%%%%%%%%Create random test signal matrix%%%%%%%%%%%%%%
    U = ones(n1,r)/sqrt(n1);
    U0 = U;
    U(mu+1:n1) = -U(mu+1:n1);
    V = ones(n2,r)/sqrt(n2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig = 1; % iid case
    indexs = 0;
    for s = ss
        indexs = indexs + 1;
        % Sdiag = linspace(6*sqrt(n1),4*sqrt(n1), r);  % tall
        % Sdiag = linspace(7*sqrt(n2),5*sqrt(n2), r);  % wide
        % Sdiag = linspace(3*sig*sqrt(n1)*log(n1)^(2/3),1.2*sig*sqrt(n1)*log(n1)^(2/3), r);  % weak signal
        % Sdiag = linspace(3*n1^(4/5),n1^(4/5), r);  % strong signal
        % Sdiag = [70,60,50];
        % Sdiag(1:r-1) = 50*n1;
        % Sdiag(r-1) = 1.1*Sdiag(r);
        Sdiag = s;
        sr = Sdiag(r);
        S = diag(Sdiag);
        M = U*S*V';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf("\n mu = %.d (%.d / %.d), s = %.d / %.d", mu, indexmu, length(mus), indexs, length(ss))
        fprintf("\n____Progress Bar____\n")
        powerfrob = zeros(reps,1);
        powerour = powerfrob;
        for k = 1:reps
            Tfrobs = zeros(m,1);
            tstats = zeros(m,1);
            for j = 1:m
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                E = normrnd(0, sig, [n1,n2]); % iid case
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mhat = M + E;
                [Uhat,Shat,Vhat] = svds(Mhat,r,'largest','MaxIterations', 1000);  
                Shat = diag(Shat);
                Shat = sqrt((Shat.^2 - (1+c)*sig^2*n2 + sqrt(((1+c)*sig^2*n2-Shat.^2).^2 - 4*c*sig^4*n2^2))/2);
                Shat = diag(Shat);
                %%%%%%%%%%%% Frobenius norm test %%%%%%%%%%%%%%
                Tfrob1 = norm(Uhat*Uhat'-U0*U0',"fro")^2+norm(Vhat*Vhat'-V*V',"fro")^2;
                Tfrobs(j) = (Tfrob1 - 2*(n1+n2-2*r)*norm(inv(Shat),"fro")^2)/(sqrt(8*(n1+n2-2*r))*norm(inv(Shat)^2,"fro"));
                %%%%%%%%%%%% Our test %%%%%%%%%%%%%%%%%%%
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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            powerfrob(k) = length(Tfrobs(Tfrobs>=qnorm))/m;
            powerour(k) = length(tstats(tstats>=qgum))/m;
            fprintf("1111")
        end    
        powerfrobmat(length(ss)*(indexmu-1) + indexs,:) = [mu, log10(s), mean(powerfrob)];
        powerourmat(length(ss)*(indexmu-1) + indexs,:) = [mu, log10(s), mean(powerour)];
    end
end
toc

save("pwr_comparison.mat","powerourmat","powerfrobmat")


