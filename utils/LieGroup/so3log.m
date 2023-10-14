function [w] = so3log(R)
%SO3LOG Log map for so(3)

theta = acos((trace(R)-1)/2);
if abs(theta) < 1e-2
    lnR = (1/2 + theta^2/12) * (R - R');
else
    lnR = theta / (2*sin(theta)) * (R - R');
end
w = [lnR(3,2); lnR(1,3); lnR(2,1)];

end

