function [shift] = phasecorr_1D(y1, y2, method, window_fn, ignore_higher_freqs)
%PHASECORR_1D
    if nargin < 3
        method = 'plane_fit';
    end    
    if nargin < 4
        window_fn = @box_window_1D;
    end
    if nargin < 5
        ignore_higher_freqs = false;
    end
    
    assert(numel(y1) == numel(y2));
    Y1 = dft_1D(y1, window_fn);
    Y2 = dft_1D(y2, window_fn);
    shift = phasecorr_dft_1D(Y1, Y2, method, ignore_higher_freqs);
end

