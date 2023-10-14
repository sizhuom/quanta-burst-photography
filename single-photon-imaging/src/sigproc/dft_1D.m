function [dft] = dft_1D(y, window_fn)
%DFT_1D
    if nargin < 2
        window_fn = @box_window_1D;
    end
    
    if ~isfloat(y)
        y = double(y);
    end
    
    if ~isequal(window_fn, @box_window_1D)
        L = numel(y);
        window = window_fn(L);
        m = mean(y);
        y = m + (y - m) .* window;
    end
    
    dft = fft(y);
end