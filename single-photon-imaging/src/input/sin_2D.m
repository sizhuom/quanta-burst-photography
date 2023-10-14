function [img] = sin_2D(L, shift, brightness, contrast, num_cycles, angle)
%GENERATE_SIN_2D
%   H, W:       height and width of resulting image
%   sx, sy:     how much to shift the image center from the original middle
%               point
%   brightness: minimum level of light
%   contrast:   0 means all-gray, 1 means black-and-white
%   num_cycles: number of cycles of the wave in the image
%   angle:      orientation of the line (should be between -180 and 180)
    if nargin < 2
        shift = [0 0];
    end
    
    if nargin < 3
        brightness = 0;
    end
    
    if nargin < 4
        contrast = 1;
    end
    
    if nargin < 5
        num_cycles = [1 1];
    end
    
    if nargin < 6
        angle = 0;
    end
    
    H = L(1); W = L(2);
    sx = shift(1); sy = shift(2);
    [X, Y] = meshgrid(-floor(W/2):-floor(W/2)+W-1, -floor(H/2):-floor(H/2)+H-1);
    
    tf = affine2d([cosd(angle) -sind(angle) 0; sind(angle) cosd(angle) 0; 0 0 1]);
    [X, Y] = transformPointsForward(tf, X - sx, Y - sy);
    
    L_diag = norm(L);
    d1 = W * abs(secd(angle)); d2 = H * abs(secd(angle));
    if max([d1 d2]) > L_diag
        d1 = H * abs(cscd(angle)); d2 = W * abs(cscd(angle));
    end
    freq = [1/d1 1/d2] .* num_cycles;
    img = normalize_minmax(sin(2 * pi * freq(1) * X) .* sin(2 * pi * freq(2) * Y));
    % at this point img is in range [0, 1]
    % map this to [0, contrast] first (minimum brightness),
    % then shift the range using the given brightness
    img = (img * contrast) + brightness;
    % no saturation implemented so far, because this is supposed to be just
    % radiance.
end