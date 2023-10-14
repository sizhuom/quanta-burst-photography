function [ Knew ] = resizeIntrinsics( K, scale )
%RESIZEINTRINSICS Modify the intrinsic matrix due to image rescaling

a = scale;
b = (1-a)/2;
Kadd = [
    a 0 b;
    0 a b;
    0 0 1;
    ];

Knew = Kadd * K;

end

