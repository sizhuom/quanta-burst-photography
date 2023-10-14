function [ flowrgb, wheelrgb ] = drawFlowHSV( flow, maxMag )
%DRAWFLOWHSV Draw flow HSV plot as images
% H: angle from x axis, counter-clockwise, 0-360
% S: magnitude, normalized by maximum flow
% V: occlusion, not implemented yet

% HSV representation of the x-y components
flowx = flow(:,:,1);
flowy = flow(:,:,2);
angle = atan2(flowy, flowx);
angle(angle<0) = angle(angle<0) + 2*pi;
h = angle / (2*pi);

mag = sqrt(flowx.^2+flowy.^2);
if nargin < 2
    maxMag = max(mag(:));
end
s = mag / maxMag;

hsv = cat(3,h,s,double(~isnan(h)));
flowrgb = hsv2rgb(hsv);

% HSV color wheel for x/y flow
if nargout > 1
    resolution = 1000;
    step = maxMag * 2 / resolution;
    [x, y] = meshgrid(-maxMag:step:maxMag);
    [theta, rho] = cart2pol(x,y);
    theta(theta<0) = theta(theta<0) + 2*pi;
    h = theta / (2*pi);
    s = min(1, rho / maxMag);
    v = rho < maxMag;
    wheelrgb = hsv2rgb(cat(3,h,s,v));
end



end

