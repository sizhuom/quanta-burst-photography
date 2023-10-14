function [yout] = align_and_average_2D(Y, shifts, avg_fn)
%ALIGN_AND_AVERAGE_1D
    [Ly, Lx, T] = size(Y);
    if nargin < 2
        shifts = zeros(T, 1);
    end
    
    if nargin < 3
        avg_fn = @mean;
    end
    
    assert(size(shifts, 1) == T && size(shifts, 2) == 2);
    Y_align = zeros(Ly, Lx, T);
    shifts = shifts - shifts(1, :);
    for t = 1:T
        Y_align(:,:,t) = translate_2D(Y(:,:,t), -shifts(t, :));
    end
    
    yout = avg_fn(Y_align, 3);
end
