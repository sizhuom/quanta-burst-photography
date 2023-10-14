function [shift] = optflow_1D(I1, I2, sigma, avg_type, Ix_threshold)
%OPTFLOW_1D
    assert(isvector(I1) && isvector(I2));
    L = numel(I1);
    
    if nargin < 3
        sigma = 0;
    end
    
    if nargin < 4
        avg_type = 'mean';
    end
    
    if nargin < 5
        Ix_threshold = 1 / L;
    end
    
    if ~isfloat(I1)
        I1 = double(I1);
    end
    if ~isfloat(I2)
        I2 = double(I2);
    end
    
    if sigma ~= 0
        g = fspecial('gaussian', [ceil(6*sigma + 1) 1], sigma);
        I1 = conv(I1, g, 'same');
        I2 = conv(I2, g, 'same');
    end    
    
    Ix = diff(I1, 1, 1);
    It = I2(1:end-1) - I1(1:end-1);
    valid_idxs = abs(Ix) >= Ix_threshold;
    It = It(valid_idxs);
    Ix = Ix(valid_idxs);
    
    if numel(Ix) > 0
        if strcmpi(avg_type, 'median')
            shift = median(-It ./ Ix);
        elseif strcmpi(avg_type, 'mean')
            shift = mean(-It ./ Ix);
        elseif strcmpi(avg_type, 'wmean')
            w = hamming(L);
            w = w(valid_idxs);
            shift = sum(w .* (-It ./ Ix)) / sum(w, 'all');
        end
    else
        shift = 0;
    end
end
