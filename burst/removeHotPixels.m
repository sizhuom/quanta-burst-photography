function im = removeHotPixels(im, dcr, thresh)
%REMOVEHOTPIXELS Remove hot pixels from binary images by interpolating from
%surrounding pixels

% Determine the hot pixel mask
hpMask = dcr > thresh;

% Process each frame
if islogical(im)
    i = mod(randi([0 7]) + 5, 9) + 1;
    h = zeros(3);
    h(i) = 1;
    imf = imfilter(im, h);
    im(hpMask) = imf(hpMask);
else
    h = ones(5) / 24;
    h(3,3) = 0;
    imf = imfilter(im, h);
    im(hpMask) = imf(hpMask);
end

