function [num_patches] = num_patches_1D(L, patch_size, patch_stride)
%NUM_PATCHES_1D
    assert(L >= patch_size);
    num_patches = 1 + floor((L - patch_size) / patch_stride);
end