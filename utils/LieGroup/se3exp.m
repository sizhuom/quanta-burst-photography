function [R, T] = se3exp(uw)
%SE3EXP Exponential map for se(3)

u = uw(1:3);
w = uw(4:6);
theta = norm(w);
if abs(theta) < 1e-2
    A = 1 - theta^2 / 6;
    B = 1/2 - theta^2 / 24;
    C = 1/6 - theta^2 / 120;
else
    A = sin(theta)/theta;
    B = (1-cos(theta))/theta^2;
    C = (1-A)/theta^2;
end
wx = skewSym(w);
wx2 = wx * wx;
R = eye(3) + A*wx + B*wx2;
V = eye(3) + B*wx + C*wx2;
T = V*u;

if nargout < 2 % Output a single matrix
    R = [R T];
end

end

