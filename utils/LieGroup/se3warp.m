function [x1, y1, Juw] = se3warp(x, y, d, RT, K)
%SE3WARP SE3 warp function
%Input:
%   x, y, d: image coordinates and inverse depth, row vector
%   RT: SE(3), 3x4 matrix
%   K: 3x3 intrinsic matrix, assuming K(3,1:2)=0 & K(3,3)=1
%Output:
%   x1, y1: warped coordinates, row vector
%   Juw: Jacobian wrt. uw, 12xN matrix

R = RT(1:3, 1:3);
T = RT(1:3, 4);

Kinv = inv(K);
p0 = Kinv * [x; y; ones(size(x))] ./ d;
p1 = R*p0 + T;
proj = [p1(1,:)./p1(3,:); p1(2,:)./p1(3,:)];
x1 = K(1,1)*proj(1,:) + K(1,2)*proj(2,:) + K(1,3);
y1 = K(2,1)*proj(1,:) + K(2,2)*proj(2,:) + K(2,3);

p13inv = 1./p1(3,:);
p13inv2 = p13inv .^ 2;
row1 = [p13inv; zeros(size(x)); -p1(1,:).*p13inv2; -p1(1,:).*p1(2,:).*p13inv2; 1+p1(1,:).^2.*p13inv2; -p1(2,:).*p13inv];
row2 = [zeros(size(x)); p13inv; -p1(2,:).*p13inv2; -1-p1(2,:).^2.*p13inv2; p1(1,:).*p1(2,:).*p13inv2; p1(1,:).*p13inv];
Juw = [K(1,1)*row1+K(1,2)*row2; K(2,1)*row1+K(2,2)*row2];

end

