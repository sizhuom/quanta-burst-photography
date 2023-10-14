function [lshifts] = get_shifts_at_lag(shifts_mat, lag)
%GET_SHIFTS_AT_LAG
    assert(ndims(shifts_mat) == 3);
    [M, N, R] = size(shifts_mat);
    assert(M == N);
    assert(abs(lag) < M);
    if lag < 0
        lag = -lag;
        % shifts_mat has to be real, so no need to conjugate
        shifts_mat = permute(shifts_mat, [2 1 3]);
    end
    
    lshifts = zeros(M - lag, R);
    for i = 1:M-lag
        lshifts(i,:) = squeeze(shifts_mat(i, i + lag, :));
    end
end

