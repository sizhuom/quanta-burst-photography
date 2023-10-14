function [ qw ] = SE32qw( SE3 )
%SE32QW Convert an SE3 matrix to qw representation
% where q is the translation vector, w is the axis-angle representation

q = SE3(1:3, 4);

R = SE3(1:3, 1:3);
theta = acos((trace(R)-1)/2);
if abs(theta) < 1e-2
    lnR = (1/2 + theta^2/12) * (R - R');
else
    lnR = theta / (2*sin(theta)) * (R - R');
end
w = [lnR(3,2); lnR(1,3); lnR(2,1)];

qw = [q; w];

end
