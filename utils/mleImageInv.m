function [S] = mleImageInv(Lambda, T, scale)
%MLEIMAGEINV Convert an linear intensity image back to sum image

if nargin < 2
    T = 1;
end
if nargin < 3
    scale = 1;
end

S = (1 - exp(-Lambda./reshape(scale,1,1,[]))) .* T;

end

