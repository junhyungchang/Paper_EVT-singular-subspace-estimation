function [U1, d] = generate_U1(U, t)
    % Input:
    %   U - n x r matrix with orthonormal columns (U' * U = I_r)
    %   d - desired 2-to-infinity Procrustes distance
    % Output:
    %   U1 - n x r matrix with orthonormal columns, delocalized,
    %        and at Procrustes distance d from U

    [n, r] = size(U);
    
    % Check that U has orthonormal columns
    if norm(U' * U - eye(r), 'fro') > 1e-10
        error('Input matrix U must have orthonormal columns.');
    end
    % Check that t is in [0,1]
    if t>1 || t<0
        error('t must be in [0,1]')
    end
    
    % Ensure that U is delocalized
    max_row_norm = max(sqrt(sum(U.^2, 2)));
    if (max_row_norm^2 * n / r) > 10  % Adjust the threshold as needed
        error('Input matrix U is not sufficiently delocalized.');
    end
    delocheck = 0;
    while delocheck == 0
        % Step 2: Create an orthonormal basis for the orthogonal complement of U
        % Perform QR decomposition on a random matrix orthogonal to U
        U_perp = null(U');
        W = U_perp * normrnd(0,1,[n - r, r]);
        [W, ~] = qr(W, 0);  % Orthonormal columns, W' * W = I_r
        W = W - U * (U' * W);  % Ensure orthogonality to U
        [W, ~] = qr(W, 0);
    
        % Step 3: Combine U and W to create U1
    %     theta = acos(1 - (d^2 * r) / n);  % Calculate rotation angle
    %     U1 = U * Q * cos(theta) + W * sin(theta);
        U1 = (1-t)*U + t*W;
    
        % Step 4: Ensure U1 has orthonormal columns
        [U1, ~] = qr(U1, 0);
    
        % Step 5: Check that U1 is delocalized
        if (tinorm(U1)^2 * n / r) < 6 % Adjust the threshold as needed
            delocheck = 1;
        end
    end
    % Step 6: Compute the Procrustes distance
    H = U' * U1;
    [X, ~, Y] = svd(H);
    R = X * Y';
    d = tinorm(U * R - U1);
%     if abs(procrustes_distance - d) > 1e-6  % Tolerance for the distance
%         % Adjust U1 to achieve the desired distance
%         U1 = U1 * (d / procrustes_distance);
%         % Re-orthonormalize U1
%         [U1, ~] = qr(U1, 0);
%     end
end
