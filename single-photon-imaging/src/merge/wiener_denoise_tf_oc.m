function [filtered_stack] = wiener_denoise_tf_oc(patch_stack, c0, patchOffsetX, patchOffsetY, window_fn)
%WIENER_DENOISE_TF_OC Denoise aligned patch stack using Wiener filtering along time axis
% Sizhuo: just filter each patch but do not merge, offset corrected
    % assert(ndims(patch_stack) == 3);
    if nargin < 5
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
    
    M = size(patch_stack, 2);
    N = size(patch_stack, 1);
    [K, L] = meshgrid(0:size(patch_stack,2)-1, 0:size(patch_stack,1)-1);
    for i = 2:D
        FT_base_shifted = FT_base .* exp(-1j*2*pi*(K*patchOffsetX(i)/M+L*patchOffsetY(i)/N));
        tmp = idft_2D(FT_base_shifted, window_fn);
        FT_frame = dft_2D(patch_stack(:,:,i), window_fn);
        Dz2 = abs(FT_frame - FT_base_shifted).^2;
        Az = Dz2 ./ (Dz2 + c * noise_variance);
        FT_merged = Az .* FT_base_shifted + (1 - Az) .* FT_frame;
        filtered_stack(:,:,i) = idft_2D(FT_merged, window_fn);
    end
end