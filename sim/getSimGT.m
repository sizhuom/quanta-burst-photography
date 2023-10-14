function [depth] = getSimGT(im,camParam)
%GETSIMGT Get ground truth depth from simulation data

im = rgb2gray(im2double(im));
depth = (1-im)*camParam.maxDepth;

end

