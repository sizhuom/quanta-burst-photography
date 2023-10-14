function [oidxs] = outlier_idxs(D, dim, sigma)
%FIND_OUTLIER_IDXS
    if nargin < 3
        sigma = 0;
    end
    
    if nargin < 2
        dim = 1;
    end
    
    if sigma ~= 0
        Ds = imgaussfilt(D, sigma);
        D = D - Ds;
    end
    
    mD = nanmedian(D, dim);
    mAD = nanmedian(abs(mD));
    oidxs = abs(mD) >= 2 * mAD;
end
