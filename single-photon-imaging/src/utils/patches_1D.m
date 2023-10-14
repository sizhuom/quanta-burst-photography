function [patches] = patches_1D(obs, patch_size, patch_stride)
%PATCHES_1D
    if nargin < 3
        patch_stride = patch_size;
    end
    
    L = size(obs, 1);
    T = size(obs, 2);
    R = size(obs, 3);
    P = num_patches_1D(L, patch_size, patch_stride);
    patches = zeros(patch_size, P, T, R);
    
    for r = 1:R
        for p = 1:P
            pstart = 1 + (p-1)*patch_stride;
            pfinish = pstart + patch_size - 1;
            patches(:,p,:,r) = obs(pstart:pfinish,:,r);
        end
    end
end