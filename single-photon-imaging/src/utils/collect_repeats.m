function [xc, yc] = collect_repeats(x, y, collect_fn)
%COLLECT_REPEATS
    if nargin < 3
        collect_fn = @nanmedian;
    end
    
    assert(numel(x) == numel(y));
    N = numel(x);
    
    [xs, is] = sort(x);
    ys = y(is);
    xc = [];
    yc = [];
    i = 1;
    while i < N
        xi = xs(i);
        yi = [ys(i)];
        for j = i+1:N
            if xs(j) ~= xi
                break;
            else
                yi = [yi ys(j)];
            end
        end
        xc = [xc xi];
        yc = [yc collect_fn(yi)];
        i = j;
%         conlog(num2str(xc(end)));
    end
end

