function R = procrustes_sol_F(U1, U2)
% % Compute the Frobenius norm optimal orthogonal Procrustes solution.
    H  = U1' * U2;
    [X,~,Y] = svd(H);
    R = X*Y';
end