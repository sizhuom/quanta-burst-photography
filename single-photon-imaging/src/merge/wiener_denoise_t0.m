function [merged] = wiener_denoise_t(patch_stack, window_fn)
%WIENER_DENOISE_T Denoise aligned patch stack using Wiener filtering along time axis
    % assert(ndims(patch_stack) == 3);
    if nargin < 2
        window_fn = @box_window_2D;
    end
    
    D = size(patch_stack, 3);
    base_frame = patch_stack(:,:,1);
    noise_variance = std(base_frame(:)) ^ 2;
    FT_base = dft_2D(base_frame, window_fn);
    FT_merged = FT_base;
    c = size(patch_stack, 1) * size(patch_stack, 2) * 2 * 8;
    for i = 2:D
        FT_frame = dft_2D(patch_stack(:,:,i), window_fn);
        Dz2 = abs(FT_frame - FT_base).^2;
        Az = Dz2 ./ (Dz2 + c * noise_variance);
        FT_merged = FT_merged + Az .* FT_base + (1 - Az) .* FT_frame;
    end
    FT_merged = FT_merged / D;
    merged = idft_2D(FT_merged, window_fn);
end