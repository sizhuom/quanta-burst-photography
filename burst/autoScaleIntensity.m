function [imScaled, imgScale] = autoScaleIntensity(im, percentile)
%AUTOSCALEINTENSITY Automatically scale the linear intensity image
%generated from real sequences

dataMat = reshape(im, [], size(im,3));
prcValue = max(prctile(dataMat, percentile));
imgScale = 1 / prcValue;
imScaled = im * imgScale;

end

