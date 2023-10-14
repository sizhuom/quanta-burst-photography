function [uw] = se3log(RT)
%SE3LOG Log map for se(3)

R = RT(1:3, 1:3);
T = RT(1:3, 4);

theta = acos((trace(R)-1)/2);
if abs(theta) < 1e-2
    lnR = (1/2 + theta^2/12) * (R - R');
else
    lnR = theta / (2*sin(theta)) * (R - R');
end
w = [lnR(3,2); lnR(1,3); lnR(2,1)];

Vinv = eye(3) - lnR/2 + (1-theta/2/tan(theta/2))/theta^2*lnR^2;
q = Vinv * T;

uw = [q; w];

end

