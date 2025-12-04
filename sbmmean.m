function M = sbmmean(B,s)
% Creates mean matrix of SBM
% B is a k by k probability matrix
% s is a vector of length k.
% k blocks, each of size s(i) by s(i).
    % Check if s is a vector of length k
    k = length(B);
    if length(s) ~= k
        error('Length of s must be equal to k.');
    end
    
    % Initialize the block matrix
    M = zeros(sum(s));
    % Generate diagonal blocks and off diagonal blocks
    start_idx = 1;
    for i = 1:k
        end_idx = start_idx + s(i) - 1;
        M(start_idx:end_idx, start_idx:end_idx) = B(i,i)*ones(s(i),s(i));
        Boff = B(i, setdiff(1:end, i));
        soff = s(setdiff(1:end, i));
        Poff = zeros(s(i), sum(s)-s(i));
        off_start_idx = 1;
        for j = 1:k-1
            off_end_idx = off_start_idx + soff(j) -1;
            Poff(:, off_start_idx:off_end_idx) = Boff(j)*ones(s(i), soff(j));
            off_start_idx = off_end_idx +1;
        end
        M(start_idx:end_idx, setdiff(1:end, start_idx:end_idx)) = Poff;
        start_idx = end_idx + 1;
    end  
end