function [depth] = getSimGTRGB(im,camParam)
%GETSIMGTRGB Get ground truth depth from simulation data

im = im2double(im);
depth = (im(:,:,2)+im(:,:,3)/128)*camParam.maxDepth;

end

