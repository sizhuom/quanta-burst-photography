clearvars; close all;
%% Generate input radiances
fprintf('Generate input radiances\n');
input_type = 'sin';

num_shots = 256;
L = 256;

initial_shift = 128;
camera_velocity = 0.01;
brightness = 0.2;
contrast = 0.8;
direction = 1;
frequency = 8/L;
smoothness = 1/3;
powerlaw_p = 1;

radiance = zeros(num_shots, L);
for s = 1:num_shots
    shift = initial_shift - camera_velocity * (s - 1);
    
    if strcmpi(input_type, 'step')
        radiance(s,:) = generate_step_1D(L, shift, brightness, contrast, direction);
    elseif strcmpi(input_type, 'sin')
        radiance(s,:) = generate_sin_1D(L, shift, brightness, contrast, frequency, direction);
    elseif strcmpi(input_type, 'smooth')
        radiance(s,:) = generate_smooth_1D(L, shift, brightness, contrast, smoothness);
    elseif strcmpi(input_type, 'power_law_fft')
        radiance(s,:) = generate_power_law_fft_1D(L, shift, brightness, contrast, powerlaw_p);
    end
end

%% Acquire binary snapshots
fprintf('Acquire binary snapshots\n');
num_trials = 256;
shot_exposure_time = 1;
b = zeros(num_trials, num_shots, L);
for s = 1:num_shots
    for t = 1:num_trials
        b(t,s,:) = get_binary_snapshot(radiance(s,:), shot_exposure_time);
    end
end

%% Preprocess binary snapshots (optional)
% fprintf('Preprocess binary snapshots\n');
% window_fn = @hamming;
% padding = L/4;
% for s = 1:num_shots
%     for t = 1:num_trials
%         obs = squeeze(b(t,s,:));
%         if padding > 0
%             obs = padarray(obs, padding, 'symmetric');
%         end
%         obs_fourier = get_1D_DFT(obs, window_fn);
%         m = abs(obs_fourier);
%         p = angle(obs_fourier);
%         m = L * log_transform_mle(m / L, L, L * shot_exposure_time);
%         obs_rescaled = get_1D_IDFT(m .* exp(1i * p), window_fn);
%         if padding > 0
%             obs_rescaled = obs_rescaled(padding+1 : padding+L);
%         end
%         b(t,s,:) = obs_rescaled;
%     end
% end

%% Extract MLE block-wise
fprintf('Extract MLE block-wise\n');
window_fn = @raised_cosine_window;
patch_size = 32;
patch_stride = 16;
num_patches = 1 + (L - patch_size) / patch_stride;
patch_weight_window = window_fn(patch_size);
mle_block_size = 32;
mle_block_stride = 16;
num_mles = 1 + ((num_shots - mle_block_size) / mle_block_stride);
y_mle = zeros(num_trials, num_mles, L);

for s = 1:num_mles
    block_start = 1 + (s-1)*mle_block_stride;
    block_end = block_start + mle_block_size - 1;
    for t = 1:num_trials
        obs_block = squeeze(b(t, block_start:block_end, :));
        obs_patches = get_patch_stack_1D(obs_block, patch_size, patch_stride);
        assert(size(obs_patches, 1) == num_patches);
        for p = 1:size(obs_patches, 1)
            patch_start = 1 + (p-1)*patch_stride;
            patch_end = patch_start + patch_size - 1;
            % TODO: proper merge procedure
            patch_mle = log_transform_mle(mean(squeeze(obs_patches(p,:,:)), 1), mle_block_size, mle_block_size * shot_exposure_time);
            y_mle(t,s,patch_start:patch_end) = y_mle(t,s,patch_start:patch_end) + reshape(patch_weight_window' .* patch_mle, [1 1 patch_size]);
        end
    end
end

%% Find pairwise shifts between MLEs
fprintf('Find pairwise shifts between MLEs\n');
window_fn = @raised_cosine_window;
shift_estimation_method = 'plane_fit';
shifts_y_mle = zeros(num_trials, num_mles, num_patches);
for s = 2:num_mles
    for t = 1:num_trials
        for p = 1:num_patches
            patch_start = 1 + (p-1)*patch_stride;
            patch_end = patch_start + patch_size - 1;
            shifts_y_mle(t, s, p) = get_phasecorr_shift_1D(y_mle(t, s-1, patch_start:patch_end), y_mle(t, s, patch_start:patch_end), window_fn, shift_estimation_method);
        end
    end
end

%% Find pairwise shifts between binary shots directly
fprintf('Find pairwise shifts between binary shots directly\n');
window_fn = @raised_cosine_window;
shift_estimation_method = 'plane_fit';
shifts_b = zeros(num_trials, num_shots, num_patches);
for s = 2:num_shots
    for t = 1:num_trials
        for p = 1:num_patches
            patch_start = 1 + (p-1)*patch_stride;
            patch_end = patch_start + patch_size - 1;
            shifts_b(t, s, p) = get_phasecorr_shift_1D(b(t, s-1, patch_start:patch_end), b(t, s, patch_start:patch_end), window_fn, shift_estimation_method);
        end
    end
end

%% Use shifts between MLEs to fix shifts between observed binary images.
% fprintf('Use shifts between MLEs to fix shifts between observed binary images\n');
% shifts_b = zeros(num_trials, num_shots, num_patches);
% ind_start = 1 + floor(mle_block_size / 2);
% ind_end = num_mles*mle_block_stride + floor(mle_block_size / 2);
% shifts_b(:, ind_start:mle_block_stride:ind_end, :) = shifts_y_mle;
% for t = 1:num_trials
%     for p = 1:num_patches
%         shifts_b(t,:,p) = smooth(shifts_b(t,:,p), 1 + 2*floor(mle_block_stride/2));
%     end
% end

true_shifts_y_mle = zeros(num_trials, num_mles, num_patches);
for m = 2:num_mles
    true_shifts_y_mle(:,m,:) = -camera_velocity * mle_block_stride;
end
mean_shifts_y_mle_centred = squeeze(mean(shifts_y_mle - true_shifts_y_mle, 1));
max_shifts_y_mle_centred = squeeze(max(shifts_y_mle - true_shifts_y_mle, [], 1));
min_shifts_y_mle_centred = squeeze(min(shifts_y_mle - true_shifts_y_mle, [], 1));
mean_shifts_y_mle_centred = squeeze(mean(mean_shifts_y_mle_centred, 2));
max_shifts_y_mle_centred = squeeze(max(max_shifts_y_mle_centred, [], 2));
min_shifts_y_mle_centred = squeeze(min(min_shifts_y_mle_centred, [], 2));

figure;
plot(mean_shifts_y_mle_centred, 'LineWidth', 2);
hold on; plot(max_shifts_y_mle_centred, '--', 'LineWidth', 2);
hold on; plot(min_shifts_y_mle_centred, '--', 'LineWidth', 2);
grid on;

true_shifts_b = zeros(num_trials, num_shots, num_patches);
for s = 2:num_shots
    true_shifts_b(:,s,:) = -camera_velocity;
end
mean_shifts_b_centred = squeeze(mean(shifts_b - true_shifts_b, 1));
max_shifts_b_centred = squeeze(max(shifts_b - true_shifts_b, [], 1));
min_shifts_b_centred = squeeze(min(shifts_b - true_shifts_b, [], 1));
mean_shifts_b_centred = squeeze(mean(mean_shifts_b_centred, 2));
max_shifts_b_centred = squeeze(max(max_shifts_b_centred, [], 2));
min_shifts_b_centred = squeeze(min(min_shifts_b_centred, [], 2));
figure;
plot(mean_shifts_b_centred, 'LineWidth', 2);
hold on; plot(max_shifts_b_centred, '--', 'LineWidth', 2);
hold on; plot(min_shifts_b_centred, '--', 'LineWidth', 2);
grid on;

% fprintf('Cheating!\n');
% shifts_b = true_shifts_b;
% fprintf('But adding some noise...\n');
% added_variance = 1e0;
% shifts_b = shifts_b + randn(size(shifts_b)) * sqrt(added_variance);

%% Generate aligned binary shots and merge
% fprintf('Aligning binary shots and merging\n');
% window_fn = @raised_cosine_window;
% b_aligned = zeros(size(b));
% % y_est_b = zeros(num_trials, L);
% patch_weight_window = window_fn(patch_size);
% for t = 1:num_trials
%     for p = 1:num_patches
%         patch_start = 1 + (p-1)*patch_stride;
%         patch_end = patch_start + patch_size - 1;
%         shift = 0;
%         for s = 1:num_shots
%             shift = shift - shifts_b(t,s,p);
%             b_aligned(t,s,patch_start:patch_end) = b_aligned(t,s,patch_start:patch_end) ...
%                                                 + reshape(patch_weight_window .* translate_1D(b(t,s,patch_start:patch_end),...
%                                                                                                 shift, 'interp1')',...
%                                                           [1 1 patch_size]);
%         end
%     end
% end
% % Merge: just the mean for now. TODO: add wiener-merge
% y_est_b = log_transform_mle(squeeze(mean(b_aligned, 2)), num_shots, num_shots * shot_exposure_time);
% mean_y_est_b = mean(y_est_b, 1);
% max_y_est_b = max(y_est_b, [], 1);
% min_y_est_b = min(y_est_b, [], 1);
% 
% % 
% t = 1:L;
% plot(radiance(1,:));
% hold on; plot(t, mean_y_est_b, 'LineWidth', 3);
% hold on; plot(t, max_y_est_b, '--', 'LineWidth', 3);
% hold on; plot(t, min_y_est_b, '--', 'LineWidth', 3);
% legend('True radiance map', 'Mean estimate after alignment', 'Max estimate', 'Min estimate');
% plot(t, radiance(1,:), t, mean(y_est_b, 1));
% 
% Y = zeros(size(radiance));
% Y_mag = zeros(size(radiance));
% Y_phase = zeros(size(radiance));
% 
% b_rescaled = zeros(size(b));
% B = zeros(size(b));
% B_mag = zeros(size(b));
% B_phase = zeros(size(b));
% 
% num_mles = 1 + ((num_shots - mle_block_size) / mle_block_stride);
% y_mle = zeros(num_trials, num_mles, L);
% Y_mle = zeros(size(y_mle));
% Y_mle_mag = zeros(size(y_mle));
% Y_mle_phase = zeros(size(y_mle));
% 
% shifts_b = zeros(num_trials, num_shots);
% shifts_y_mle = zeros(num_trials, num_mles);
% window_fn = @box_window;
% shift_estimation_method = 'plane_fit';
% 
% for s = 1:num_shots
%     
%     Y(s,:) = fft(radiance(s,:));
%     Y_mag(s,:) = abs(Y(s,:));
%     Y_phase(s,:) = unwrap(angle(Y(s,:)));
%     
%     for t = 1:num_trials
%         b(t,s,:) = get_binary_snapshot(radiance(s,:), shot_exposure_time);
%         B(t,s,:) = get_1D_DFT(b(t,s,:), window_fn);
%         B_mag(t,s,:) = abs(B(t,s,:));
%         B_phase(t,s,:) = unwrap(angle(B(t,s,:)));
%         b_rescaled(t,s,:) = get_1D_IDFT(L * log_transform_mle(B_mag(t,s,:) / L, L, L)  .* exp(1i * B_phase(t,s,:)), window_fn);
%     end
%     
%     if s >= mle_block_size && mod(s - mle_block_size, mle_block_stride) == 0
%         i_block = 1 + ((s - mle_block_size) / mle_block_stride);
%         b_block = b(:, s-mle_block_size+1:s, :);
%         for tt = 1:num_trials
%             y_mle(tt, i_block, :) = log_transform_mle(mean(b_block(tt,:,:), 2), mle_block_size, mle_block_size);
%             y_mle(tt, i_block, :) = mean(b_block(tt,:,:), 2);
%             Y_mle(tt, i_block, :) = fft(y_mle(tt, i_block, :));
%             Y_mle_mag(tt, i_block, :) = abs(Y_mle(tt, i_block, :));
%             Y_mle_phase(tt, i_block, :) = unwrap(angle(Y_mle(tt, i_block, :)));
%         end
%     end
% end
% 
% for s = 2:num_shots
%     for t = 1:num_trials
% %         shifts_b(t, s) = get_phasecorr_shift_1D(b(t, s-1, :), b(t, s, :), window_fn, shift_estimation_method);
%         shifts_b(t, s) = get_phasecorr_shift_1D(b_rescaled(t, s-1, :), b_rescaled(t, s, :), window_fn, shift_estimation_method);
%     end
% end
% 
% for s = 2:num_mles
%     for t = 1:num_trials
%         shifts_y_mle(t, s) = get_phasecorr_shift_1D(y_mle(t, s-1, :), y_mle(t, s, :), window_fn, shift_estimation_method);
%     end
% end
% 
% true_shift_direct = camera_velocity;
% true_shift_mle = camera_velocity * mle_block_stride;
% mean_shift_direct = mean(shifts_b, 'all');
% mean_shift_mle = mean(shifts_y_mle, 'all');
% rmse_direct = sqrt(mean((shifts_b - true_shift_direct) .^ 2, 'all'));
% rmse_mle = sqrt(mean((shifts_y_mle - true_shift_mle) .^ 2, 'all'));
% disp([true_shift_direct mean_shift_direct rmse_direct]);
% disp([true_shift_mle mean_shift_mle rmse_mle]);
% 
% y_reconstructed_block_mle = zeros(num_trials, L);
% for t = 1:num_trials
%     y_reconstructed_block_mle(t,:) = y_mle(t,1,:);
%     shift = 0;
%     for s = 2:num_mles
%         shift = shift - shifts_y_mle(t,s);
%         y_reconstructed_block_mle(t,:) = y_reconstructed_block_mle(t,:) + translate_1D(squeeze(y_mle(t,s,:)), shift, 'interp1');
%     end
% end
% y_reconstructed_block_mle = y_reconstructed_block_mle / num_mles;
% 
% y_reconstructed_pairwise = zeros(num_trials, L);
% for t = 1:num_trials
%     y_reconstructed_pairwise(t,:) = b_rescaled(t,1,:);
%     shift = 0;
%     for s = 2:num_shots
%         shift = shift - shifts_b(t,s);
%         y_reconstructed_pairwise(t,:) = y_reconstructed_pairwise(t,:) + translate_1D(squeeze(b_rescaled(t,s,:)), shift, 'interp1');
%     end
% end
% y_reconstructed_pairwise = y_reconstructed_pairwise / num_shots;

% figure;
% 
% subplot(1, 2, 1);
% plot(y);
% title('Original radiance map');
% 
% subplot(1, 2, 2);
% plot(squeeze(mean(b, [1 2])));
% title('Mean photon detections by sensor array');
% 
% y_mle = zeros(num_trials, L);
% Y_mle = zeros(size(y_mle));
% Y_mle_mag = zeros(size(y_mle));
% Y_mle_phase = zeros(size(y_mle));
% for t = 1:num_trials
%     y_mle(t,:) = get_flux_mle(mean(b(t,:,:), 2), mle_block_size, mle_block_size);
%     Y_mle(t,:) = fft(y_mle(t,:));
%     Y_mle_mag(t,:) = abs(Y_mle(t,:));
%     Y_mle_phase(t,:) = unwrap(angle(Y_mle(t,:)));
% end
% 
% mean_magY = mean(Y_mle_mag, 1);
% std_magY = std(Y_mle_mag, 0, 1);
% snr_magY = 20 * log10(mean_magY ./ std_magY);
% 
% mean_phaseY = mean(Y_mle_phase, 1);
% std_phaseY = std(Y_mle_phase, 0, 1);
% snr_phaseY = 20 * log10(abs(mean_phaseY) ./ std_phaseY);
% 
% mean_magB = squeeze(mean(B_mag, [1 2]));
% std_magB = squeeze(std(B_mag, 0, [1 2]));
% snr_magB = 20 * log10(mean_magB ./ std_magB);
% 
% mean_phaseB = squeeze(mean(B_phase, [1 2]));
% std_phaseB = squeeze(std(B_phase, 0, [1 2]));
% snr_phaseB = 20 * log10(abs(mean_phaseB) ./ std_phaseB);
% 
% figure;
% 
% subplot(2, 4, 1);
% plot(y);
% title('Spatial domain signal');
% 
% subplot(2, 4, 2);
% plot(mean_magY);
% xlabel('(Scaled) Frequency');
% ylabel('Mean MLE Coefficient Magnitude');
% title('MLE Magnitude spectrum');
% 
% subplot(2, 4, 3);
% plot(snr_magY);
% xlabel('(Scaled) Frequency');
% ylabel('SNR (dB) of MLE magnitude spectrum');
% title('SNR vs frequency (magnitude)');
% 
% subplot(2, 4, 4);
% plot(snr_phaseY);
% xlabel('(Scaled) Frequency');
% ylabel('SNR (dB) of MLE phase spectrum');
% title('SNR vs frequency (phase)');
% 
% subplot(2, 4, 7);
% plot(snr_magB);
% xlabel('(Scaled) Frequency');
% ylabel('SNR (dB) of direct magnitude spectrum');
% title('SNR vs frequency (magnitude)');
% 
% subplot(2, 4, 8);
% plot(snr_phaseB);
% xlabel('(Scaled) Frequency');
% ylabel('SNR (dB) of direct phase spectrum');
% title('SNR vs frequency (phase)');
% 
% subplot(2, 4, 5);
% scatter(Y_mag, snr_magB, 'DisplayName', 'Raw binary');
% hold on; scatter(Y_mag, snr_magY, 'filled', 'DisplayName', 'MLE');
% xlabel('Magnitude of Fourier coefficient');
% ylabel('SNR (dB) of magnitude spectrum');
% legend('Location', 'southeast');
% title('SNR vs coefficient magnitude (mag)');
% 
% subplot(2, 4, 6);
% scatter(Y_mag, snr_phaseB, 'DisplayName', 'Raw binary');
% hold on; scatter(Y_mag, snr_phaseY, 'filled', 'DisplayName', 'MLE');
% xlabel('Magnitude of Fourier coefficient');
% ylabel('SNR (dB) of phase spectrum');
% legend('Location', 'southeast');
% title('SNR vs coefficient magnitude (phase)');

% Plot distribution of coefficients.
% num_hist_bins = 32;
% hist_magB = zeros(num_hist_bins, L);
% hist_phaseB = zeros(num_hist_bins, L);
% for i = 1:L
%     hist_magB(:, i) = histcounts(B_mag(:,:,i), num_hist_bins);
%     hist_phaseB(:,i) = histcounts(B_phase(:,:,i), num_hist_bins);
% end

% figure;
% 
% subplot(1, 3, 1);
% plot(Y_mag);
% xlabel('(scaled) frequency');
% ylabel('coefficient magnitude');
% title('Magnitude spectrum of original radiance map');
% 
% subplot(1, 3, 2);
% surf(hist_magB);
% xlabel('(scaled) frequency');
% ylabel('(scaled) coefficient magnitude');
% zlabel('(scaled) count');
% title('Magnitude spectrum distribution');
% 
% subplot(1, 3, 3);
% surf(hist_phaseB);
% xlabel('(scaled) frequency');
% ylabel('(scaled) coefficient magnitude');
% zlabel('(scaled) count');
% title('Phase spectrum distribution');
% suptitle('Fourier coefficient distributions');