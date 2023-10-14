function [y] = idft_1D(Y, window_fn)
%GET_1D_IDFT
    y = real(ifft(Y));
    if ~isequal(window_fn, @box_window_1D)
        L = numel(Y);
        window = window_fn(L);
        y = y(:) ./ window(:);
    end
end