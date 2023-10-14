function [yout] = align_and_average_1D(Y, shifts, avg_fn)
%ALIGN_AND_AVERAGE_1D
    [L, T] = size(Y);
    if nargin < 2
        shifts = zeros(T, 1);
    end
    
    if nargin < 3
        avg_fn = @mean;
    end
    
    assert(isvector(shifts) && numel(shifts) == T);
    Y_align = zeros(L, T);
    shifts = shifts - shifts(1);
    for t = 1:T
        Y_align(:,t) = translate_1D(Y(:,t), -shifts(t));
    end
    
    yout = avg_fn(Y_align, 2);
end
