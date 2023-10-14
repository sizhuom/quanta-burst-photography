function [y] = step_1D(L, shift, brightness, contrast, direction)
%STEP_1D Step function with various parameters.
    if nargin < 2
        shift = L / 2;
    end
    
    if nargin < 3
        brightness = 0;
    end
    
    if nargin < 4
        contrast = 1;
    end
    
    if nargin < 5
        direction = 1;
    end
    
    x = 0:L-1;
    y = brightness + zeros(L, 1);
    if direction > 0
        y(x >= shift) = y(x >= shift) + contrast;
    elseif direction < 0
        y(x < shift) = y(x < shift) + contrast;
    end
end