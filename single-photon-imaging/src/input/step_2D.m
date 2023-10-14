function [img] = generate_step_2D(H, W, sx, sy, brightness, contrast, angle)
%GENERATE_STEP_2D
%   H, W:       height and width of resulting image
%   sx, sy:     how much to shift the image center from the original middle
%               point
%   brightness: minimum level of light
%   contrast:   0 means all-gray, 1 means black-and-white
%   angle:      angle of the line (should be between -180 and 180)
    img = zeros(H, W);
    [X, Y] = meshgrid(-floor(W/2):-floor(W/2)+W-1, -floor(H/2):-floor(H/2)+H-1);
    if abs(angle - 90) < 0.1
        img(X < 0)  = -contrast/2;
        img(X >= 0) =  contrast/2;
    elseif abs(angle + 90) < 0.1
        img(X < 0)  =  contrast/2;
        img(X >= 0) = -contrast/2;
    else
        m = tan(deg2rad(angle));
        img(m * (X - sx) <=  (Y - sy)) = -contrast/2;
        img(m * (X - sx) > (Y - sy)) =  contrast/2;
    end
    % at this point img is in range [-contrast/2, contrast/2]
    % map this to [0, contrast] first (minimum brightness),
    % then shift the range using the given brightness
    img = img + contrast/2 + brightness;
    % no saturation implemented so far, because this is supposed to be just
    % radiance.
end