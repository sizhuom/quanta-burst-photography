function [flux_mle] = reconstruct_flux_from_binary_snapshots_mle(snapshots, exposure_time, merge_method, align_parameters, merge_parameters)
%RECONSTRUCT_FLUX_FROM_BINARY_SNAPSHOTS_MLE Max likelihood estimate of
%incoming flux from binary snapshots assuming aligned shots
    assert(ndims(snapshots) == 3);  % assume a stack of 2d images.
    num_shots = size(snapshots, 3);
    if strcmp(merge_method, 'average')
        s_mean = mean(snapshots, 3);
    elseif strcmp(merge_method, 'wiener')
        original_shots = snapshots;
        if align_parameters.use_denoised_shots_for_alignment
            snapshots = double(snapshots);
            for i = 1:size(snapshots, 3)
                snapshots(:,:,i) = wiener_filter_spatial(snapshots(:,:,i), 8, merge_parameters.spatial_wienerfilt_freq_tuning_exponent);
            end
        end
        image_patches = get_patch_stack(snapshots, align_parameters);
        if ~align_parameters.cheat
            patch_shifts = fourier_align(image_patches, align_parameters);
        else
            [num_patches, ~, ~, D] = size(image_patches);
            shifts = get_linear_translations(align_parameters.true_linear_velocity, D);
            patch_shifts = repmat(reshape(shifts, [1 D 2]), num_patches, 1, 1);
        end
        
        image_patches = get_patch_stack(original_shots, align_parameters);
        if max(abs(patch_shifts(:))) ~= 0
            aligned_patch_stack = apply_patch_shifts(double(image_patches), -patch_shifts, align_parameters);
        else
            aligned_patch_stack = image_patches;
        end
        s_mean = wiener_merge(aligned_patch_stack, merge_parameters);
    end
    flux_mle = get_flux_mle(s_mean, num_shots, exposure_time);
end