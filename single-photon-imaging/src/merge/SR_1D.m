function [y, xq] = SR_1D(Y, shifts, zoom, repeats_collect_fn)
%SR_1D
    if nargin < 3
        zoom = 1;
    end
    
    if nargin < 4
        repeats_collect_fn = @nanmedian;
    end
    
    [L, T] = size(Y);
    assert(isvector(shifts) && numel(shifts) == T);
    s = shifts(:, ones(L, 1))';
    xs = (1:L)';
    xs = xs(:, ones(T, 1));
    xs = xs - s; % coordinate so far = true + shift, therefore true = so far - shift
    xf = reshape(xs, [], 1);
    yf = reshape(Y, [], 1);
    [xf, yf] = collect_repeats(xf, yf, repeats_collect_fn);
    x0 = min(xf, [], 'all') + (1/zoom);
    x1 = max(xf, [], 'all') - (1/zoom);
    xq = x0:(1/zoom):x1;
    y = interp1(xf, yf, xq, 'spline');
end