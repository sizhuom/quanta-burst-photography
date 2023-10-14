function [Mlr] = lowrank(M, k, zeronan)
%LOWRANK
    if nargin < 3
        zeronan = true;
    end
    
    if zeronan
        M(isnan(M)) = 0;
    end
    
    [U, S, V] = svd(M);
    Mlr = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
end

