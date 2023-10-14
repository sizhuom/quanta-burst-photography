function [y] = plmag_1D(L, shift, brightness, contrast, p, randomness, random_seed)
%plmag_1D
    % shift only works exactly for integer values; others will cause
    % artifacts.
    if nargin < 2
        shift = 0;
    end
    
    if nargin < 3
        brightness = 0;
    end
    
    if nargin < 4
        contrast = 1;
    end
    
    if nargin < 5
        p = 1;
    end
    
    if nargin < 6
        randomness = 'uniform';
    end
    
    if nargin < 7
        random_seed = 'default';
    end
        
    f = fftshift((0:L-1)' - L/2) * (1 / L) * 2 * pi;
    Y_mag = 1.0 ./ ((abs(f) + 1).^p);

    incoming_rng = rng;
    rng(random_seed);
    if strcmpi(randomness, 'uniform')
        Y_phase = 2 * pi * rand(size(Y_mag));
    elseif strcmpi(randomness, 'gaussian') || strcmpi(randomness, 'normal')
        Y_phase = wrapTo2Pi(2 * pi * randn(size(Y_mag)));
    end
    rng(incoming_rng);

    Y_phase = Y_phase - f*shift;
    Y = Y_mag .* exp(1i * Y_phase);
    y = idft_1D(Y, @box_window_1D);
    y = normalize_minmax(y);
    y = brightness + contrast * y;
end