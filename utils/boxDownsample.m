function imd = boxDownsample(im, boxSize, normalize)
%BOXDOWNSAMPLE Downsample by taking average of each nxn patch
% 210208: minor fix, add option to not normalize
if nargin < 3
    normalize = true;
end
scaleFilter = ones(boxSize);
if normalize
    scaleFilter = scaleFilter / boxSize^2;
end
scaleCenter = floor((boxSize+1)/2);
imd = imfilter(im, scaleFilter);
imd = imd(scaleCenter:boxSize:end,scaleCenter:boxSize:end,:);

end

