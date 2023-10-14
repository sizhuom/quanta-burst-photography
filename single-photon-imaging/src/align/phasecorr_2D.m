function [shift] = phasecorr_2D(y1, y2, method, window_fn)
%PHASECORR_2D
    if nargin < 3
        method = 'plane_fit';
    end
    if nargin < 4
        window_fn = @box_window_2D;
    end
    
    [M1, N1] = size(y1);
    [M2, N2] = size(y2);
    assert(M1 == M2 && N1 == N2);
    
    Y1 = dft_2D(y1, window_fn);
    Y2 = dft_2D(y2, window_fn);
    shift = phasecorr_dft_2D(Y1, Y2, method);
end
