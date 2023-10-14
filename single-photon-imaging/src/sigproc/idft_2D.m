function [y] = idft_2D(Y, window_fn)
%IDFT_2D
    y = real(ifft2(Y));
    if ~isequal(window_fn, @box_window_2D)
        [H, W] = size(Y);
        window = window_fn(H, W);
        y = y ./ window;
    end
end

