function [bias_shifts_y_mle, iqr_shifts_y_mle, bias_shifts_b_from_ymle, iqr_shifts_b_from_ymle] ...
                                = simulate_1D_reconstruction(input_type, num_shots, signal_length, num_trials, ...
                                                             camera_velocity, brightness, contrast, shot_exposure_time, ...
                                                             patch_size, patch_stride, window_fn, shift_estimation_method, ...
                                                             mle_block_size, mle_block_stride)
%SIMULATE_1D_RECONSTRUCTION
    %% Generate input radiances
    fprintf('Generate input radiances\n');
    L = signal_length;

    initial_shift = 0;
    if strcmpi(input_type, 'step')
        initial_shift = L / 2;
    end
    direction = 1;
    frequency = 8/L;
    smoothness = 1/3;
    powerlaw_p = 1;

    radiance = zeros(L, num_shots);
    for s = 1:num_shots
        shift = initial_shift - camera_velocity * (s - 1);

        if strcmpi(input_type, 'step')
            radiance(:,s) = generate_step_1D(L, shift, brightness, contrast, direction);
        elseif strcmpi(input_type, 'sin')
            radiance(:,s) = generate_sin_1D(L, shift, brightness, contrast, frequency, direction);
        elseif strcmpi(input_type, 'smooth')
            radiance(:,s) = generate_smooth_1D(L, shift, brightness, contrast, smoothness);
        elseif strcmpi(input_type, 'power_law_fft')
            radiance(:,s) = generate_power_law_fft_1D(L, shift, brightness, contrast, powerlaw_p);
        end
    end
    
    %% Acquire binary snapshots
    fprintf('Acquire binary snapshots\n');
    binary_obs = zeros(L, num_shots, num_trials);
    for s = 1:num_shots
        for t = 1:num_trials
            binary_obs(:,s,t) = get_binary_snapshot(radiance(:,s), shot_exposure_time);
        end
    end

    %% Preprocess binary snapshots (optional)
%     fprintf('Preprocess binary snapshots\n');
%     padding = L/4;
%     for s = 1:num_shots
%         for t = 1:num_trials
%             obs = binary_obs(:,s,t);
%             if padding > 0
%                 obs = padarray(obs, padding, 'symmetric');
%             end
%             obs_fourier = get_1D_DFT(obs, window_fn);
%             m = abs(obs_fourier);
%             p = angle(obs_fourier);
%             m = L * log_transform_mle(m / L, L, L * shot_exposure_time);
%             obs_rescaled = get_1D_IDFT(m .* exp(1i * p), window_fn);
%             if padding > 0
%                 obs_rescaled = obs_rescaled(padding+1 : padding+L);
%             end
%             binary_obs(:,s,t) = obs_rescaled;
%         end
%     end

    %% Set up patch size etc. for later steps (common for many)
    num_patches = get_num_patches_1D(L, patch_size, patch_stride);

    %% Find pairwise shifts between binary shots directly
    % fprintf('Find pairwise shifts between binary shots directly\n');
    % shifts_pairwise_b_direct = zeros(num_patches, num_shots, num_trials);
    % for p = 1:num_patches
    %     patch_start = 1 + (p-1)*patch_stride;
    %     patch_end = patch_start + patch_size - 1;
    %     for s = 2:num_shots
    %         for t = 1:num_trials
    %             prev_patch = binary_obs(patch_start:patch_end,s-1,t);
    %             this_patch = binary_obs(patch_start:patch_end, s, t);
    %             shifts_pairwise_b_direct(p,s,t) = get_phasecorr_shift_1D(prev_patch, this_patch, window_fn, shift_estimation_method);
    %         end
    %     end
    % end

    %% Extract MLE block-wise
    fprintf('Extract MLE block-wise\n');
    patch_weight = window_fn(patch_size);
    num_mles = 1 + ((num_shots - mle_block_size) / mle_block_stride);
    y_mle = zeros(L, num_mles, num_trials);
    for s = 1:num_mles
        block_start = 1 + (s-1)*mle_block_stride;
        block_end = block_start + mle_block_size - 1;
        for t = 1:num_trials
            obs_block = binary_obs(:, block_start:block_end, t);
            obs_patches = get_patch_stack_1D(obs_block, patch_size, patch_stride);
            assert(size(obs_patches, 3) == num_patches);
            for p = 1:num_patches
                patch_start = 1 + (p-1)*patch_stride;
                patch_end = patch_start + patch_size - 1;
                % TODO: proper merge procedure
                patch_mle = log_transform_mle(mean(obs_patches(:,:,p), 2), mle_block_size, mle_block_size * shot_exposure_time);
                y_mle(patch_start:patch_end,s,t) = y_mle(patch_start:patch_end,s,t) + patch_weight .* patch_mle;
            end
        end
    end

    %% Find pairwise shifts between MLEs
    fprintf('Find pairwise shifts between MLEs\n');
    shifts_pairwise_y_mle = zeros(num_patches, num_mles, num_trials);
    for p = 1:num_patches
        patch_start = 1 + (p-1)*patch_stride;
        patch_end = patch_start + patch_size - 1;
        for s = 2:num_mles
            for t = 1:num_trials
                prev_mle = y_mle(patch_start:patch_end, s-1, t);
                this_mle = y_mle(patch_start:patch_end, s, t);
                shifts_pairwise_y_mle(p,s,t) = get_phasecorr_shift_1D(prev_mle, this_mle, window_fn, shift_estimation_method);
            end
        end
    end

    %% Use shifts between MLEs to fix shifts between observed binary images.
    fprintf('Use shifts between MLEs to fix shifts between observed binary images\n');
    shifts_pairwise_b_from_ymle = zeros(num_patches, num_shots, num_trials);
    for m = 1:num_mles
        mleshift = shifts_pairwise_y_mle(:,m,:);
        ind_end = mle_block_size + (m-1)*mle_block_stride;
        ind_start = ind_end - mle_block_size + 1;
        distributed_shift = (1.0 / mle_block_stride) * repmat(mleshift, 1, mle_block_size, 1);    
        if m < 3
            retention_factor = 0;
        else
            retention_factor = (mle_block_size - mle_block_stride) / mle_block_size;
        end        
        shifts_pairwise_b_from_ymle(:, ind_start:ind_end, :) = retention_factor * shifts_pairwise_b_from_ymle(:, ind_start:ind_end, :) ...
                                                                + (1 - retention_factor) * distributed_shift;
    end

    %% Compare obtained shifts to the true values.
    true_shifts_pairwise_y_mle = zeros(num_patches, num_mles);
    for m = 2:num_mles
        true_shifts_pairwise_y_mle(:,m) = -camera_velocity * mle_block_stride;
    end
    true_shifts_pairwise_b = zeros(num_patches, num_shots);
    for s = 2:num_shots
        true_shifts_pairwise_b(:,s) = -camera_velocity;
    end

    t_mle = mle_block_size:mle_block_stride:L;
    shifts_y_mle_collected = reshape(permute(shifts_pairwise_y_mle, [2 1 3]), [num_mles num_patches*num_trials]);
    % plot_shift_distribution(t_mle, shifts_y_mle_collected, true_shifts_pairwise_y_mle(1,:)', ...
    %                         sprintf('Blockwise MLE(%d,%d) shift', mle_block_size, mle_block_stride), ...
    %                         'sample index', 'shift (px)');
    bias_shifts_y_mle = nanmedian(shifts_y_mle_collected - repmat(true_shifts_pairwise_y_mle(1,:)', 1, num_patches*num_trials), 'all');
    iqr_shifts_y_mle = median(get_interquantile_range(shifts_y_mle_collected, 0.05, 0.95), 'all');

    t_direct = 1:L;
    % shifts_b_direct_collected = reshape(permute(shifts_pairwise_b_direct, [2 1 3]), [num_shots num_patches*num_trials]);
    % plot_shift_distribution(t_direct, shifts_b_direct_collected, true_shifts_pairwise_b(1,:)', ...
    %                         'Directly recovered shift', 'sample index', 'shift (px)');
    % bias_shifts_b_direct = nanmedian(shifts_b_direct_collected - repmat(true_shifts_pairwise_b(1,:)', 1, num_patches*num_trials), 'all');
    % iqr_shifts_b_direct = median(get_interquantile_range(shifts_b_direct_collected, 0.05, 0.95), 'all');

    shifts_b_from_ymle_collected = reshape(permute(shifts_pairwise_b_from_ymle, [2 1 3]), [num_shots num_patches*num_trials]);
    % plot_shift_distribution(t_direct, shifts_b_from_ymle_collected, true_shifts_pairwise_b(1,:)', ...
    %                         sprintf('Evenly distributed MLE(%d,%d)-shift', mle_block_size, mle_block_stride), ...
    %                         'sample index', 'shift (px)');
    bias_shifts_b_from_ymle = nanmedian(shifts_b_from_ymle_collected - repmat(true_shifts_pairwise_b(1,:)', 1, num_patches*num_trials), 'all');
    iqr_shifts_b_from_ymle = median(get_interquantile_range(shifts_b_from_ymle_collected, 0.05, 0.95), 'all');    
end