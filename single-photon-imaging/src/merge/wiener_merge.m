function [merged] = wiener_merge(aligned_patch_stack, merge_parameters)
%WIENER_MERGE Merge image stack using Wiener filtering
%   Assume everything is already aligned.
    W = merge_parameters.merged_size_x;
    H = merge_parameters.merged_size_y;
    window_fn = merge_parameters.window_fn;
    use_whole_stack = merge_parameters.use_whole_stack;

    [~, patch_size_y, patch_size_x, D] = size(aligned_patch_stack);
    assert(mod(H, patch_size_y/2) == 0);
    assert(mod(W, patch_size_x/2) == 0);
    
    merged = zeros(H, W);
    % Weight makes sure we only take good parts from each reconstruction.
    % For the cosine window, we've also ensured the weights sum to 1 when
    % added up in a sliding window fashion.
    patch_weight = get_2D_window(window_fn, patch_size_y, patch_size_x);
    i_patch = 0; % Arrrrrr!!
    for i = 1:patch_size_y/2:H-patch_size_y+1
        for j = 1:patch_size_x/2:W-patch_size_x+1
            i_patch = i_patch + 1;
            aligned_patches = squeeze(aligned_patch_stack(i_patch,:,:,:));
            if use_whole_stack
                merged_patch = wiener_filter_depth(aligned_patches, window_fn);
            else
                merged_patch = squeeze(aligned_patches(:,:,1));
            end
            merged(i:i+patch_size_y-1, j:j+patch_size_x-1) = merged(i:i+patch_size_y-1, j:j+patch_size_x-1) + (patch_weight .* merged_patch);
        end
    end
    
    if merge_parameters.do_spatial_postfilt
        merged = wiener_filter_spatial(merged, 8 / D, merge_parameters.spatial_wienerfilt_freq_tuning_exponent);
    end
end