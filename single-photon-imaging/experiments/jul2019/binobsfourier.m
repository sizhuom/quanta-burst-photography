clearvars;
close all;

dft_window_fn = @box_window_1D;

L = 256;
input_fn = @sin_1D;
offset = 0;
if isequal(input_fn, @step_1D)
    offset = L/2;
end

brightness = 0.2;
contrast = 0.8;

camera_velocity = -0.01 * 5;
T = 256;

%% Generate input signal.
Y = zeros(L, T);
Y_F = zeros(size(Y));
for t = 1:T
    shift = offset - camera_velocity*(t-1);
    Y(:,t) = input_fn(L, shift, brightness, contrast);
%     Y_F(:,t) = dft_1D(Y(:,t), dft_window_fn);
end

%% Generate SPAD samples from input signal.
R = 32;
B = false(L, T, R);
B_F = zeros(size(B));
for r = 1:R
    for t = 1:T
        B(:,t,r) = get_spad_shot(Y(:,t));
%         B_F(:,t,r) = dft_1D(B(:,t,r), dft_window_fn);
    end
end

%% Generate blockwise MLE, assuming no motion.
mle_bsize = 32;
mle_bstride = 16;
M = num_patches_1D(T, mle_bsize, mle_bstride);
Y_mle = zeros(L, M, R);
Y_mle_F = zeros(size(Y_mle));
Y_bavg = zeros(L, M);
Y_bavg_F = zeros(size(Y_bavg));
for m = 1:M
    bstart = 1 + (m-1)*mle_bstride;
    bend = bstart + mle_bsize - 1;
    for r = 1:R
        b_block = B(:,bstart:bend,r);
        Y_mle(:,m,r) = radiance_mle_spad(mean(b_block, 2), mle_bsize, mle_bsize);
%         Y_mle_F(:,m,r) = dft_1D(Y_mle(:,m,r), dft_window_fn);
    end
    Y_bavg(:,m) = mean(Y(:,bstart:bend), 2);
%     Y_bavg_F(:,m) = dft_1D(Y_bavg(:,m), dft_window_fn);
end

%% Get magnitude and phase spectra separately.
% Only applies to sin_1D.
% if isequal(input_fn, @sin_1D)
%     % Remove frequencies not of interest.
%     Y_F(3:end-1,:) = nan;
%     B_F(3:end-1,:,:) = nan;
%     Y_mle_F(3:end-1,:,:) = nan;
%     Y_bavg_F(3:end-1,:,:) = nan;
% end
% 
% % Select one frequency and one pair of instances.
% Y_F = Y_F([2], 1:2);
% B_F = B_F([2], 1:2, :);
% Y_mle_F = Y_mle_F([2], 1:2, :);
% Y_bavg_F = Y_bavg_F([2], 1:2, :);
% 
% Y_Fmag = abs(Y_F);
% Y_Fph = angle(Y_F);
% 
% B_Fmag = abs(B_F);
% B_Fph = angle(B_F);
% 
% Y_mle_Fmag = abs(Y_mle_F);
% Y_mle_Fph = angle(Y_mle_F);
% Y_bavg_Fmag = abs(Y_bavg_F);
% Y_bavg_Fph = angle(Y_bavg_F);

% Y_Fmag = repmat(Y_Fmag, 1, 1, R);
% Y_Fph = repmat(Y_Fph, 1, 1, R);
% Y_bavg_Fmag = repmat(Y_bavg_Fmag, 1, 1, R);
% Y_bavg_Fph = repmat(Y_bavg_Fph, 1, 1, R);
% Y_Fmag = reshape(Y_Fmag, [], 1);
% B_Fmag = reshape(B_Fmag, [], 1);
% Y_mle_Fmag = reshape(Y_mle_Fmag, [], 1);
% Y_bavg_Fmag = reshape(Y_bavg_Fmag, [], 1);
% Y_Fph = reshape(Y_Fph, [], 1);
% B_Fph = reshape(B_Fph, [], 1);
% Y_mle_Fph = reshape(Y_mle_Fph, [], 1);
% Y_bavg_Fph = reshape(Y_bavg_Fph, [], 1);

% figure;
% subplot(1, 2, 1);
% scatter(Y_Fmag, B_Fmag);
% hold on; scatter(Y_bavg_Fmag, Y_mle_Fmag); hold off;
% subplot(1, 2, 2);
% scatter(Y_Fph, B_Fph);
% hold on; scatter(Y_bavg_Fph, Y_mle_Fph); hold off;

%% Get shifts between all pairs of observations/MLEs
true_shifts_B = zeros(T, T, R);
shifts_B = zeros(T, T, R);
for r = 1:R
    for t1 = 1:T
        for t2 = 1:T
            shifts_B(t1,t2,r) = phasecorr_1D(B(:,t1,r), B(:,t2,r), 'inv_ft', dft_window_fn);
            true_shifts_B(t1,t2,r) = -camera_velocity * (t2 - t1);
        end
    end
end

true_shifts_Ymle = zeros(M, M, R);
shifts_Ymle = zeros(M, M, R);
for r = 1:R
    for m1 = 1:M
        for m2 = 1:M
            shifts_Ymle(m1,m2,r) = phasecorr_1D(Y_mle(:,m1,r), Y_mle(:,m2,r), 'inv_ft', dft_window_fn);
            true_shifts_Ymle(m1,m2,r) = -camera_velocity * (m2 - m1) * mle_bstride;
        end
    end
end

%% Clean up shift estimates by enforcing the vector sum identity.
shifts_B_smooth = zeros(size(shifts_B));
shifts_Ymle_smooth = zeros(size(shifts_Ymle));
for r = 1:R
    shifts_B_smooth(:,:,r) = force_vectorsum_equality(shifts_B(:,:,r), 1);
    shifts_Ymle_smooth(:,:,r) = force_vectorsum_equality(shifts_Ymle(:,:,r), 1);
end


%% Find variance in shift estimates as a function of lag.
%  lag goes from 1 to T/M - 1.
shifts_B_markers = zeros(T-1, 3);
true_shifts_B_markers = zeros(T-1, 1);
for lag = 1:T-1
    lshifts = get_shifts_at_lag(shifts_B, lag);
    shifts_B_markers(lag, 1) = min(lshifts, [], 'all');
    shifts_B_markers(lag, 2) = median(lshifts, 'all');
    shifts_B_markers(lag, 3) = max(lshifts, [], 'all');
    true_lshifts = get_shifts_at_lag(true_shifts_B, lag);
    true_shifts_B_markers(lag) = true_lshifts(1);
end

shifts_B_smooth_markers = zeros(T-1, 3);
true_shifts_B_markers = zeros(T-1, 1);
for lag = 1:T-1
    lshifts = get_shifts_at_lag(shifts_B_smooth, lag);
    shifts_B_smooth_markers(lag, 1) = min(lshifts, [], 'all');
    shifts_B_smooth_markers(lag, 2) = median(lshifts, 'all');
    shifts_B_smooth_markers(lag, 3) = max(lshifts, [], 'all');
    true_lshifts = get_shifts_at_lag(true_shifts_B, lag);
    true_shifts_B_markers(lag) = true_lshifts(1);
end

shifts_Ymle_markers = zeros(M-1, 3);
true_shifts_Ymle_markers = zeros(M-1, 1);
for lag = 1:M-1
    lshifts = get_shifts_at_lag(shifts_Ymle, lag);
    shifts_Ymle_markers(lag, 1) = min(lshifts, [], 'all');
    shifts_Ymle_markers(lag, 2) = median(lshifts, 'all');
    shifts_Ymle_markers(lag, 3) = max(lshifts, [], 'all');
    true_lshifts = get_shifts_at_lag(true_shifts_Ymle, lag);
    true_shifts_Ymle_markers(lag) = true_lshifts(1);
end

shifts_Ymle_smooth_markers = zeros(M-1, 3);
true_shifts_Ymle_markers = zeros(M-1, 1);
for lag = 1:M-1
    lshifts = get_shifts_at_lag(shifts_Ymle_smooth, lag);
    shifts_Ymle_smooth_markers(lag, 1) = min(lshifts, [], 'all');
    shifts_Ymle_smooth_markers(lag, 2) = median(lshifts, 'all');
    shifts_Ymle_smooth_markers(lag, 3) = max(lshifts, [], 'all');
    true_lshifts = get_shifts_at_lag(true_shifts_Ymle, lag);
    true_shifts_Ymle_markers(lag) = true_lshifts(1);
end

figure;
subplot(2, 2, 1);
lags = 1:T-1;
plot(lags, shifts_B_markers(:,1), '--', 'LineWidth', 2, 'DisplayName', 'minimum shift');
hold on; 
plot(lags, shifts_B_markers(:,2), 'LineWidth', 3, 'DisplayName', 'median shift');
plot(lags, shifts_B_markers(:,3), '--', 'LineWidth', 2, 'DisplayName', 'maximum shift');
plot(lags, true_shifts_B_markers, '--k', 'LineWidth', 2, 'DisplayName', 'true shift');
hold off;
xlabel('Lag');
ylabel('Estimated shift (px)');
title('Shift estimates from binary shots, as function of lag');
legend;

subplot(2, 2, 2);
lags = (1:M-1)*mle_bstride;
plot(lags, shifts_Ymle_markers(:,1), '--', 'LineWidth', 2, 'DisplayName', 'minimum shift');
hold on;
plot(lags, shifts_Ymle_markers(:,2), '-o', 'LineWidth', 3, 'DisplayName', 'median shift');
plot(lags, shifts_Ymle_markers(:,3), '--', 'LineWidth', 2, 'DisplayName', 'maximum shift');
plot(lags, true_shifts_Ymle_markers, '--k', 'LineWidth', 2, 'DisplayName', 'true shift');
hold off;
xlabel('Lag');
ylabel('Estimated shift (px)');
title('Shift estimates from blockwise MLEs, as function of lag');

subplot(2, 2, 3);
lags = 1:T-1;
plot(lags, shifts_B_smooth_markers(:,1), '--', 'LineWidth', 2, 'DisplayName', 'minimum shift');
hold on; 
plot(lags, shifts_B_smooth_markers(:,2), 'LineWidth', 3, 'DisplayName', 'median shift');
plot(lags, shifts_B_smooth_markers(:,3), '--', 'LineWidth', 2, 'DisplayName', 'maximum shift');
plot(lags, true_shifts_B_markers, '--k', 'LineWidth', 2, 'DisplayName', 'true shift');
hold off;
xlabel('Lag');
ylabel('Estimated shift (px)');
title('Smoothed shift estimates from binary shots, as function of lag');

subplot(2, 2, 4);
lags = (1:M-1)*mle_bstride;
plot(lags, shifts_Ymle_smooth_markers(:,1), '--', 'LineWidth', 2, 'DisplayName', 'minimum shift');
hold on;
plot(lags, shifts_Ymle_smooth_markers(:,2), '-o', 'LineWidth', 3, 'DisplayName', 'median shift');
plot(lags, shifts_Ymle_smooth_markers(:,3), '--', 'LineWidth', 2, 'DisplayName', 'maximum shift');
plot(lags, true_shifts_Ymle_markers, '--k', 'LineWidth', 2, 'DisplayName', 'true shift');
hold off;
xlabel('Lag');
ylabel('Estimated shift (px)');
title('Smoothed shift estimates from blockwise MLEs, as function of lag');

