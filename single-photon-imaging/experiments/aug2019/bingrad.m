clearvars;
close all;

N = 10000;
R = 10;

b0 = logical(randi(1, N, N, R));
b1 = logical(randi(1, size(b0)));

conlog('Binary xor');
% gb = false(size(b0));
tic;
gb = xor(b0, b1);
% for r = 1:R
%     gb(:, :, r) = xor(b0(:, :, r), b1(:, :, r));
%     gb = xor(b0, b1);
% end
% gb = xor(b0, b1);
toc;

d0 = rand(N, N, R);
d1 = rand(N, N, R);
% gd = zeros(size(d0));
conlog('Double sub');
tic;
% for r = 1:R
%     gd(:, :, r) = d1(:, :, r) - d0(:, :, r);
%     gd = d1 - d0;
% end
gd = d1 - d0;
toc;

% 
% conlog('Binary diff')
% tic;
% % b0sy1 = circshift(b0, -1, 1);
% % dyb0 = xor(b0, b0sy1);
% [dy2b0, sdy2b0] = mag_grad_2D_b(b0, 1);
% % dyb0 = false(N-1, N);
% % for t = 1:N-1
% %     dxb0(:, t) = xor(b0(:,t), b0(:,t+1));
% %     dyb0(t, :) = xor(b0(t,:), b0(t+1,:));
% % end
% toc;
% 
% conlog('Double diff');
% tic;
% % d0sx1 = circshift(d0, -1, 2);
% % dxd0 = d0 - d0sx1;
% % dxd0 = diff(d0, 1, 2);
% dyd0 = diff(d0, 1, 1);
% toc;

% Y0 = load_grayscale_img('cameraman.tif');
% b0 = get_spad_shot(Y0);
% Y1 = translate_2D(Y0, [0.1 0.05]);
% b1 = get_spad_shot(Y1);
% 
% T = 1000;
% 
% conlog('double optflow');
% tic;
% for t = 1:T
%     s0 = optflow_2D(b0, b1);
% end
% toc;
% 
% conlog('binary optflow');
% tic;
% for t = 1:T
%     s1 = optflow_2D_b(b0, b1);
% end
% toc;