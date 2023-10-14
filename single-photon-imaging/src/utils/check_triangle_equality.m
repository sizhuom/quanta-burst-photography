function [is_satisfied] = check_triangle_equality(D, eps)
%CHECK_TRIANGLE_EQUALITY
    if nargin < 2
        eps = 1e-6;
    end
    
    [M, N] = size(D);
    assert(M == N);
    is_satisfied = true;
    for i = 1:M
        for j = 1:M
            for k = 1:M
                if abs(D(i,j) + D(j,k) - D(i,k)) > eps
                    is_satisfied = false;
                end
            end
        end
    end        
end

