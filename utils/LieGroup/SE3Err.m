function [trErr] = SE3Err(RT0, RT1)
%SE3ERR Compute the translational and rotational error from SE(3) matrices

tErr = norm(RT0(1:3,4) - RT1(1:3,4));
rErr = norm(so3log(RT0(1:3,1:3)' * RT1(1:3,1:3)));

trErr = [tErr; rErr];

end

