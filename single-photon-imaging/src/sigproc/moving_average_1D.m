function [Ys] = moving_average_1D(Y, dim, span, stride, avg_fn)
%TEMPORAL_MA
    if nargin < 2 || numel(dim) == 0
        dim = 1;
    end
    
    if nargin < 4
        stride = span;
    end
    
    if nargin < 5
        if islogical(Y)
            avg_fn = @majority;
        else
            avg_fn = @mean;
        end
    end
        
    T = size(Y, dim);
    Ts = num_patches_1D(T, span, stride);
    sz = size(Y);
    sz(dim) = Ts;
    if islogical(Y)
        Ys = false(sz);
    else
        Ys = zeros(sz);
    end
    
    all_idxs = repmat({':'}, 1, ndims(Y));
    for t = 1:Ts
        block_idxs.type = '()';
        block_idxs.subs = all_idxs;
        block_idxs.subs(dim) = {1+(t-1)*stride:(t-1)*stride+span};
        s_idx.type = '()';
        s_idx.subs = all_idxs;
        s_idx.subs(dim) = {t};
        Ys = subsasgn(Ys, s_idx, avg_fn(subsref(Y, block_idxs), dim));
    end        
end