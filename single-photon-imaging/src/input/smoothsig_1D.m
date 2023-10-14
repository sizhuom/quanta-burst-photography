function [y] = smoothsig_1D(L, shift, brightness, contrast, smoothness, randomness, random_seed)
%SMOOTHSIG_1D Generate some smooth signal in 1D.
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
        smoothness = 1/2;
    end
    
    if nargin < 6
        randomness = 'uniform';
    end
    
    if nargin < 7
        random_seed = 'default';
    end
    
    k = 50;
    incoming_rng = rng;
    rng(random_seed);
    if strcmpi(randomness, 'uniform')
        x = rand(4 * k * L, 1);
    elseif strcmpi(randomness, 'gaussian') || strcmpi(randomness, 'normal')
        x = randn(4 * k * L, 1);
    end
    xs = smooth(x, max(1, floor(smoothness*3*k*L)));
    xs = normalize_minmax(xs(k*L+1:3*k*L));
    rng(incoming_rng);
    
    t = 0:2*k*L-1;
    t_shift = mod(t - k*shift, 2*k*L);
    xr = interp1(t, xs, t_shift, 'spline', 0.5);
    i = (k/2)*L+1:k:3*(k/2)*L;
    y = xr(i);
    y = normalize_minmax(y);
    y = brightness + contrast * y;    
end