function [Ab] = diag_bandify(A, k, pval)
%DIAG_BANDIFY
    assert(size(A, 1) == size(A, 2));
    [X, Y] = meshgrid(1:size(A, 1));
    Ab = A;
    Ab(abs(X - Y) > k) = pval;
end

