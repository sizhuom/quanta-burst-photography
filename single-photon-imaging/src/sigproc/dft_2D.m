function [Y] = dft_2D(y, window_fn)
%DFT_2D
    if nargin < 2
        window_fn = @box_window_2D;
    end
    
    if ~isfloat(y)
        y = double(y);
    end
    
    if ~isequal(window_fn, @box_window_2D)
        [H, W] = size(y);
        window = window_fn(H, W);
        m = mean(y, 'all');
        y = m + (y - m) .* window;
    end
    
    Y = fft2(y);
end
