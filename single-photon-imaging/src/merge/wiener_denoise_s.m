function [filtered] = wiener_denoise_s(img, const_scale, freq_exponent)
%WIENER_DENOISE_S Wiener denoising in space
    if nargin < 2
        const_scale = 2;
    end
    
    if nargin < 3
        freq_exponent = 8;
    end
    
    [H, W] = size(img);
    noise_variance = std(img(:)) ^ 2;
    dft = fft2(img);
    mag2 = abs(dft) .^ 2;
    % Attenuate high-frequency coefficients as f^p (f is normalized to
    % [0,1], so the weight is always <= 1)
    % (this weight is in the denominator of the filter, so high value =>
    % more attenuation)
    % Higher the power, the more the weight of the lower frequencies is 
    % driven to zero, meaning that coefficient is not attenuated. 
    % This results in sharper but noisier images.
    ww = (linspace(0, 1, W/2)') .^ freq_exponent;
    wh = (linspace(0, 1, H/2)') .^ freq_exponent;
    weight = wh * ww';
    weight = [weight, fliplr(weight); ...
              flipud(weight), flipud(fliplr(weight))];
    c = H^2 * W^2 * const_scale * weight;
    % dft = dft .* (mag2 ./ (mag2 + c * noise_variance));
    deno = max(eps, mag2 + c * noise_variance);
    dft = dft .* (mag2 ./ deno); %Sizhuo: quick fix to avoid divide-by-0
    filtered = real(ifft2(dft));
    filtered(filtered < 0) = 0;
    filtered(filtered > 1) = 1;
end