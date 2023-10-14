function [filtered_stack] = wiener_denoise_tf(patch_stack, c0, window_fn)
%WIENER_DENOISE_TF Denoise aligned patch stack using Wiener filtering along time axis
% Sizhuo: just filter each patch but do not merge
    % assert(ndims(patch_stack) == 3);
    if nargin < 3
        window_fn = @box_window_2D;
    end
    
    D = size(patch_stack, 3);
    base_frame = patch_stack(:,:,1);
    %noise_variance = std(base_frame(:)) ^ 2;
    noise_variance = max(eps, std(base_frame(:)) ^ 2); %Sizhuo: quick fix to avoid divide-by-0
    FT_base = dft_2D(base_frame, window_fn);
    filtered_stack = zeros(size(patch_stack));
    filtered_stack(:,:,1) = base_frame;
    c = size(patch_stack, 1) * size(patch_stack, 2) * 2 * c0;
    for i = 2:D
        FT_frame = dft_2D(patch_stack(:,:,i), window_fn);
        Dz2 = abs(FT_frame - FT_base).^2;
        Az = Dz2 ./ (Dz2 + c * noise_variance);
        FT_merged = Az .* FT_base + (1 - Az) .* FT_frame;
        filtered_stack(:,:,i) = idft_2D(FT_merged, window_fn);
    end
end