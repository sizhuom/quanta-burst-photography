function I = normalizeImg(D, validRange)
%NORMALIZEIMG Normalize an image to the range [0, 1]

if nargin < 2
    validRange = [min(D(:)) max(D(:))];
end

I = (D - validRange(1)) / (validRange(2) - validRange(1));
I(I < 0) = 0;
I(I > 1) = 0;

end

