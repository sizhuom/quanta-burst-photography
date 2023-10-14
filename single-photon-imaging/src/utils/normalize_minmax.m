function [norm_Y] = normalize_minmax(Y)
    norm_Y = (Y - min(Y(:))) / (max(Y(:)) - min(Y(:)));
end
