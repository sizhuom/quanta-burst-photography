function [y] = sin_1D(L, shift, brightness, contrast, num_cycles)
%SIN_1D A sine wave with various parameters
    if nargin < 2
        shift = 0;
    end
    
    if nargin < 3
        brightness = 0;
    end
    
    if nargin < 4
        contrast = 1;
    end
    
    if nargin < 5
        num_cycles = 1;
    end
    
    x = (0:L-1)';
    frequency = (1/L) * num_cycles;
    y = brightness + contrast/2 + (contrast/2) * sin((2 * pi * frequency * (x - shift)));
end