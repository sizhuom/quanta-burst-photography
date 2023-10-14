clearvars;
close all;

%% Generate input radiance and camera shifts.
input_type = 'sin';
L = 256;    % length of signal (spatial extent)
offset = 0;
if strcmpi(input_type, 'step')
    offset = L/2;
end
brightness = 0.2;
contrast = 0.8;

camera_velocity = -0.01 * 5;
T = 256; % number of shots, equivalently amount of 't'ime

Y = zeros(L, T); % using 'y' for radiance
true_shifts_pf = zeros(T, 1); % pf: previous frame
for t = 1:T
    shift = offset - camera_velocity * (t-1);
    if strcmpi(input_type, 'sin')
        Y(:,t) = sin_1D(L, shift, brightness, contrast);
    elseif strcmpi(input_type, 'step')
        Y(:,t) = step_1D(L, shift, brightness, contrast);
    elseif strcmpi(input_type, 'smooth')
        Y(:,t) = smoothsig_1D(L, shift, brightness, contrast, smoothness);
    elseif strcmpi(input_type, 'powerlaw')
        Y(:,t) = plmag_1D(L, shift, brightness, contrast);
    end
    
    if t >= 2
        true_shifts_pf(t) = -camera_velocity;
    end
end

%% Set up patch size for motion/shift estimation
%  me: motion estimation, s: spatial, t: temporal, p: patch, b: block
%  patch is always spatial, block always temporal
me_psize = L;
me_pstride = L;
me_pweight = box_window_1D(me_psize);
% me_pweight = raised_cos_window(me_psize);
P_me = num_patches_1D(L, me_psize, me_pstride); % P: number of patches

%% Use true_shifts to re-align the ground truth (baseline for possible reconstruction performance)
% realign_bcond = 'constant'; % boundary condition
% realign_extrap = 'edge_vals';
% Y_realign_global = zeros(size(Y));
% Y_realign_pwise = zeros(size(Y));
% shift = 0;
% for t = 1:T
%     shift = shift + true_shifts_pf(t);
%     for p = 1:P_me
%         pstart = 1 + (p-1)*me_pstride;
%         pend = pstart + me_psize - 1;
%         if shift >= 0
%             nbrhood = Y(pstart : min(pend+ceil(shift), L), t);
%             unshifted_nbrhood = translate_1D(nbrhood, -shift, realign_bcond, realign_extrap);
%             unshifted_patch = unshifted_nbrhood(1:me_psize);
%         else
%             nbrhood = Y(max(1, pstart-ceil(shift)) : pend, t);
%             unshifted_nbrhood = translate_1D(nbrhood, -shift, realign_bcond, realign_extrap);
%             unshifted_patch = unshifted_nbrhood(end-me_psize+1:end);
%         end
% %         unshifted_patch = translate_1D(Y(pstart:pend,t), -shift, realign_bcond, realign_extrap);
%         Y_realign_pwise(pstart:pend,t) = Y_realign_pwise(pstart:pend,t) ...
%                                          + me_pweight .* unshifted_patch;
%     end
%     Y_realign_global(:,t) = translate_1D(Y(:,t), -shift, realign_bcond, realign_extrap);
% end
% 
% figure;
% plot(Y(:,1), '--', 'LineWidth', 6, 'DisplayName', 'original');
% title('True input at different instants');
% for t = round(T/7):round(T/7):T-7
%     shift = sum(true_shifts_pf((1:t)));
%     hold on; plot(Y(:,t), 'LineWidth', 2, 'DisplayName', sprintf('t = %d, shift = %.2f', t, shift));
% end
% legend;
% 
% figure;
% plot(Y(:,1), '--', 'LineWidth', 4, 'DisplayName', 'original');
% title(sprintf('Patchwise re-aligned true input assuming true shifts -- %s boundary', realign_bcond));
% for t = round(T/7):round(T/7):T-7
%     shift = sum(true_shifts_pf((1:t)));
%     hold on; plot(Y_realign_pwise(:,t), 'LineWidth', 2, 'DisplayName', sprintf('t = %d, shift = %.2f', t, shift));
% end
% legend;
% 
% figure;
% plot(Y(:,1), '--', 'LineWidth', 4, 'DisplayName', 'original');
% title(sprintf('Globally re-aligned true input assuming true shifts -- %s boundary', realign_bcond));
% for t = round(T/7):round(T/7):T-7
%     shift = sum(true_shifts_pf((1:t)));
%     hold on; plot(Y_realign_global(:,t), 'LineWidth', 2, 'DisplayName', sprintf('t = %d, shift = %.2f', t, shift));
% end
% legend;

%% Get binary observations
R = 32; % number of 'r'uns
B = false(L, T, R); % 'b'inary observations
for t = 1:T
    for r = 1:R
        B(:,t,r) = get_spad_shot(Y(:,t));
    end
end

%% Get shifts between binary frames.
% pshifts_pf_mat = zeros(P_me, T, T, R);
% dft_window_fn = @box_window;
% for p = 1:P_me
%     pstart = 1 + (p-1)*me_pstride;
%     pend = pstart + me_psize - 1;
%     for r = 1:R
%         for t1 = 1:T
%             for t2 = 1:T
%                 if t1 ~= t2
%                     p1 = B(pstart:pend, t1, r);
%                     p2 = B(pstart:pend, t2, r);
%                     pshifts_pf_mat(p, t1, t2, r) = phasecorr_1D(p1, p2, 'plane_fit', dft_window_fn, false, true);
%                 end
%             end
%         end
%     end
% end

%% Extract MLE blockwise, assuming no motion
%  s: spatial, t: temporal, p: patch, b: block
%  patch is always spatial, block always temporal
mle_psize = me_psize;
mle_pstride = me_pstride;
mle_pweight = me_pweight;
P_mle = num_patches_1D(L, mle_psize, mle_pstride);
mle_bsize = 32;
mle_bstride = 16;
M = num_patches_1D(T, mle_bsize, mle_bstride); % number of 'M'LEs
true_shifts_mle = zeros(M, 1);
Y_mle = zeros(L, M, R);
for m = 1:M
    t_start = 1 + (m-1)*mle_bstride;
    t_end = t_start + mle_bsize - 1;
    for r = 1:R
        block = B(:, t_start:t_end, r);
        block_patches = patches_1D(block, mle_psize, mle_pstride);
        for p = 1:P_mle
            pstart = 1 + (p-1)*mle_pstride;
            pend = pstart + mle_psize - 1;
            patch_mle = radiance_mle_spad(mean(block_patches(:,p,:), 3), ...
                                          mle_bsize, mle_bsize);
            Y_mle(pstart:pend,m,r) = Y_mle(pstart:pend,m,r) ...
                                     + patch_mle .* mle_pweight;
        end
    end
    
    if m >= 2
        true_shifts_mle(m) = -camera_velocity * mle_bstride;
    end
end

% Do the same for the true radiances, to get an idea about the difference
% made by the observation model.
Y_bavg = zeros(L, M);
for m = 1:M
    t_start = 1 + (m-1)*mle_bstride;
    t_end = t_start + mle_bsize - 1;
    block = Y(:, t_start:t_end);
    block_patches = patches_1D(block, mle_psize, mle_pstride);
    for p = 1:P_mle
        pstart = 1 + (p-1)*mle_pstride;
        pend = pstart + mle_psize - 1;
        patch_avg = mean(block_patches(:,p,:), 3);
        Y_bavg(pstart:pend, m) = Y_bavg(pstart:pend, m) ...
                                 + patch_avg .* mle_pweight;
    end
end

%% Extract shifts between MLEs and between the true block averages.
pshifts_mle = zeros(P_me, M, R);
pshifts_bavg = zeros(P_me, M);
dft_window_fn = @box_window_1D;
for p = 1:P_me
    pstart = 1 + (p-1)*me_pstride;
    pend = pstart + me_psize - 1;
    for m = 2:M
        for r = 1:R
            prev_mle = Y_mle(pstart:pend, m-1, r);
            this_mle = Y_mle(pstart:pend, m, r);
            pshifts_mle(p,m,r) = phasecorr_1D(prev_mle, this_mle, 'plane_fit', dft_window_fn);
        end
        prev_bavg = Y_bavg(pstart:pend, m-1);
        this_bavg = Y_bavg(pstart:pend, m);
        pshifts_bavg(p, m) = phasecorr_1D(prev_bavg, this_bavg, 'plane_fit', dft_window_fn);
    end
end

t = mle_bsize:mle_bstride:T;
shifts_bavg_collected = permute(pshifts_bavg, [2 1]);
median_shifts_bavg = nanmedian(shifts_bavg_collected, 2);
max_shifts_bavg = max(shifts_bavg_collected, [], 2);
min_shifts_bavg = min(shifts_bavg_collected, [], 2);
figure;
plot(t, median_shifts_bavg, '-o', 'LineWidth', 3, 'DisplayName', 'median shift');
hold on;
plot(t, max_shifts_bavg, '-o', 'DisplayName', 'max');
plot(t, min_shifts_bavg, '-o', 'DisplayName', 'min');
plot(t, true_shifts_mle, '--k', 'LineWidth', 3, 'DisplayName', 'True value');
title('Recovered shifts (true block average)');
legend;
hold off;

shifts_mle_collected = reshape(permute(pshifts_mle, [2 1 3]), M, P_me*R);
median_shifts_mle = nanmedian(shifts_mle_collected, 2);
max_shifts_mle = max(shifts_mle_collected, [], 2);
min_shifts_mle = min(shifts_mle_collected, [], 2);
upq_shifts_mle = quantile(shifts_mle_collected, 0.95, 2);
downq_shifts_mle = quantile(shifts_mle_collected, 0.05, 2);
figure;
plot(t, median_shifts_mle, '-o', 'LineWidth', 3, 'DisplayName', 'median shift');
hold on;
plot(t, max_shifts_mle, '-o', 'DisplayName', 'max');
plot(t, min_shifts_mle, '-o', 'DisplayName', 'min');
plot(t, upq_shifts_mle, '-o', 'DisplayName', '95 percentile');
plot(t, downq_shifts_mle, '-o', 'DisplayName', '5 percentile');
plot(t, true_shifts_mle, '--k', 'LineWidth', 3, 'DisplayName', 'True value');
title('Recovered shifts (MLE)');
legend;
hold off;

%% Convert shifts in MLEs to shifts in observations.
% pshifts_pf_mle = zeros(P_me, T, R);
% for m = 1:M
%     mleshift = pshifts_mle(:,m,:);
%     i1 = mle_bsize + (m-1)*mle_bstride;
%     i0 = i1 - mle_bsize + 1;
%     distr_shift = (1.0 / mle_bstride) * repmat(mleshift, 1, mle_bsize, 1); % distributed shift
%     if m < 3
%         rf = 0; % retention factor
%     else
%         rf = (mle_bsize - mle_bstride) / mle_bsize;
%     end
%     pshifts_pf_mle(:, i0:i1, :) = rf * pshifts_pf_mle(:, i0:i1, :) ...
%                                  + (1 - rf) * distr_shift;
% end

%% Choose which shift to apply to each observation
% pshifts_pf = pshifts_pf_mle;
% pshifts_pf = zeros(P_me, T, R);
% pshifts_pf = repmat(reshape(true_shifts_pf, [1 T 1]), P_me, 1, R);

% t = 1:T;
% shifts_pf_collected = reshape(permute(pshifts_pf, [2 1 3]), T, P_me*R);
% median_shifts_pf = median(shifts_pf_collected, 2);
% max_shifts_pf = max(shifts_pf_collected, [], 2);
% min_shifts_pf = min(shifts_pf_collected, [], 2);
% upq_shifts_pf = quantile(shifts_pf_collected, 0.95, 2);
% downq_shifts_pf = quantile(shifts_pf_collected, 0.05, 2);
% figure;
% plot(t, median_shifts_pf, 'LineWidth', 3, 'DisplayName', 'median shift');
% hold on;
% plot(t, max_shifts_pf, 'DisplayName', 'max');
% plot(t, min_shifts_pf, 'DisplayName', 'min');
% plot(t, upq_shifts_pf, 'DisplayName', '95 percentile');
% plot(t, downq_shifts_pf, 'DisplayName', '5 percentile');
% plot(t, true_shifts_pf, '--k', 'LineWidth', 3, 'DisplayName', 'True value');
% title('Recovered shifts (per frame)');
% legend;
% hold off;

%% Apply obtained shifts to get a single image.
% realign_bcond = 'constant';
% realign_extrap = 'edge_vals';
% B_patches_aligned = zeros(me_psize, P_me, T-mle_bsize+mle_bstride, R);
% for r = 1:R
%     for p = 1:P_me
%         pstart = 1 + (p-1)*me_pstride;
%         pend = pstart + me_psize - 1;
%         shift = 0;
%         t0 = mle_bsize - mle_bstride;
%         for t = t0+1:T
%             shift = shift + pshifts_pf(p, t, r);
%             if shift >= 0
%                 nbrhood = double(B(pstart : min(pend+ceil(shift), L), t, r));
%                 unshifted_nbrhood = translate_1D(nbrhood, -shift, realign_bcond, realign_extrap);
%                 unshifted_patch = unshifted_nbrhood(1:me_psize);
%             else
%                 nbrhood = double(B(max(1, pstart-ceil(-shift)) : pend, t, r));
%                 unshifted_nbrhood = translate_1D(nbrhood, -shift, realign_bcond, realign_extrap);
%                 unshifted_patch = unshifted_nbrhood(end-me_psize+1:end);
%             end
%             B_patches_aligned(:,p,t-t0,r) = unshifted_patch;
%         end
%     end
% end
% B_patches_merged = permute(mean(B_patches_aligned, 3), [1 2 4 3]);
% Y0_est = zeros(L, R);
% for r = 1:R
%     for p = 1:P_me
%         pstart = 1 + (p-1)*me_pstride;
%         pend = pstart + me_psize - 1;
%         patch_est = radiance_mle_spad(B_patches_merged(:,p,r), T, T);
%         Y0_est(pstart:pend,r) = Y0_est(pstart:pend,r) ...
%                                 + patch_est .* me_pweight;
%     end
% end
% 
% figure;
% t = 0:L-1;
% plot(repmat(t, 1, R), reshape(Y0_est, T*R, 1), '.', 'DisplayName', 'reconstructions');
% hold on; plot(t, Y(:,1), '--y', 'LineWidth', 6, 'DisplayName', 'original');
% legend;
% title('Y_{est} from re-aligned binary shots');

%% Use the shifts to localize the observations in a common reference frame (the first frame).
% locs_obs = zeros(me_psize, P_me, T, R);
% for p = 1:P_me
%     pstart = 1 + (p-1)*me_pstride;
%     pend = pstart + me_psize - 1;
%     shift = 0;
%     for t = 1:T
%         shift = shift - pshifts_pf(p, t, :);
%         locs_obs(:,p,t,:) = repmat((pstart:pend)', 1, 1, 1, R) ...
%                             + repmat(reshape(shift, [1 1 1 R]), ...
%                                      me_psize, 1, 1, 1);
%     end
% end
% 
% B_patches = patches_1D(B, me_psize, me_pstride);
% B_patches = reshape(B_patches, [], R);
% locs_obs = reshape(locs_obs, [], R);
% [locs_obs, I] = sort(locs_obs, 1);
% B_patches = B_patches(I);