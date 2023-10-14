clearvars;
close all;

L = 32;
t = (0:L-1)';

a = 0.1;
b = 0.03;
% y1 = sin_1D(L, 0, 0, 1, 0.1);
y1 = a + b * t;

delta = 5;
% y2 = sin_1D(L, delta, 0, 1, 0.1);
y2 = a + b * (t - delta);

plot(t, y1, t, y2);
legend('y1', 'y2');

d_ld = line_delta_1D(y1, y2);

% bhat1 = (y1(end) - y1(1)) / (t(end) - t(1));
% bhat2 = (y2(end) - y2(1)) / (t(end) - t(1));
% bhat = geomean([bhat1 bhat2]);
% dhat = (y1(1) - y2(1)) / bhat;
% 
dpc = phasecorr_1D(y1, y2, 'plane_fit', @box_window_1D);