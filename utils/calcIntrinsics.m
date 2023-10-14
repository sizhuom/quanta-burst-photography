function K = calcIntrinsics(resol, hfov)
%CALCINTRINSICS Calculate intrinsics from resolution and horizontal FOV
% Notice that i-th pixel expands from x=i-0.5 to x=i+0.5, and is reprensented by
% the ray at i
%INPUT:
% resol: resolution of the image [height width]
% hfov: horizontal FOV in degrees
%OUTPUT:
% K: intrinsic matrix

K = [resol(2)*0.5/tan(hfov/2/180*pi) 0 resol(2)/2+0.5;
    0 resol(2)*0.5/tan(hfov/2/180*pi) resol(1)/2+0.5;
    0 0 1];

end

