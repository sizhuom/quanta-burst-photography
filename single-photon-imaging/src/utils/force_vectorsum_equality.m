function [Ds] = force_vectorsum_equality(D, alpha)
%FORCE_VECTORSUM_EQUALITY
    if nargin < 2
        alpha = 1;
    end
    
    [M, N, P] = size(D);
    assert(M == N && P <= 2);
    Ds = D;
    for r = 1:M
        for c = r+1:M
            d_rc = zeros(M-1, P);
            i = 1;
            d_rc(i, :) = D(r, c, :);
            for k = 1:r-1
                i = i + 1;
                d_rc(i, :) = reshape(D(k, c, :) - D(k, r, :), [1 P]);
            end
            for k = r+1:c-1
                i = i + 1;
                d_rc(i, :) = reshape(D(r, k, :) + D(k, c, :), [1 P]);
            end
            for k = c+1:M
                i = i + 1;
                d_rc(i, :) = reshape(D(r, k, :) - D(c, k, :), [1 P]);
            end
            valid_idxs = (sum(isnan(d_rc), 2) == 0);
            if sum(valid_idxs) > 0
                Ds(r, c, :) = alpha * reshape(nanmedian(d_rc, 1), [1 1 P]) + (1 - alpha) * Ds(r, c, :);
                Ds(c, r, :) = -Ds(r, c, :);
            end
        end
    end
end