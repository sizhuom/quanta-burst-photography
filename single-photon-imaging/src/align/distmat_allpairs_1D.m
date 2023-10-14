function D = distmat_allpairs_1D(Y, dist_fn)
%DISTMAT_ALLPAIRS_1D
    [~, T, R, S] = size(Y);

    D = zeros(T, T, R, S);
    for t1 = 1:T
        D(t1, t1, :) = 0;
        for t2 = [1:t1-1 t1+1:T]
            for r = 1:R
                for s = 1:S
                    D(t1, t2, r, s) = dist_fn(Y(:, t1, r, s), Y(:, t2, r, s));
                end
            end
        end
    end        
end