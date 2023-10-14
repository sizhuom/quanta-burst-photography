clearvars;
% close all;

dft_window_fn = @box_window_2D;
shift_est_method = 'optflow';

L = [256 256];
brightness = 0.2;
contrast = 0.8;

vs = [0.1 0.05];
T = 256;
max_shift = 2;
half_shift_window = min(T, ceil(max_shift / norm(vs)));

%% Generate the input signal and the true shifts.
conlog('Generating the input signal');
D0 = [zeros(1, 2); cumsum(vs(ones(T-1,1), :), 1)];
[t2, t1] = meshgrid(1:T);
D0_x = D0(:,1); D0_y = D0(:,2);
D_true = zeros(T, T, 2);
D_true(:,:,1) = D0_x(t2) - D0_x(t1); D_true(:,:,2) = D0_y(t2) - D0_y(t1);

Y = zeros([L T]);
Y0 = imresize(load_grayscale_img('cameraman.tif'), L);
for t = 1:T
    Y(:,:,t) = translate_2D(Y0, D0(t,:));
%     Y(:,:,t) = sin_2D(L, D0(t,:), brightness, contrast, [1 1], 0);
end

Y_noalign = align_and_average_2D(Y, zeros(T, 2), @mean);
Y_aligned = align_and_average_2D(Y, D0, @mean);

%% Generate baseline shift matrix
% conlog('Estimating shifts from the true signal');
% tic;
% c_ofsigma = norm(vs);
% D_Y = nan(T, T, 2);
% if strcmpi(shift_est_method, 'phasecorr')
%     Y_F = zeros(size(Y));
%     for t = 1:T
%         Y_F(:,:,t) = dft_2D(Y(:,:,t), dft_window_fn);
%     end
%     for t1 = 1:T
%         Y_F_t1 = Y_F(:,:,t1);
%         D_Y(t1,t1,:) = 0;
%         for t2 = [max(1,t1-max_shift):t1-1 t1+1:min(T,t1+max_shift)]
%             Y_F_t2 = Y_F(:,:,t2);
%             D_Y(t1, t2, :) = phasecorr_dft_2D(Y_F_t1, Y_F_t2, 'plane_fit');
%         end
%     end
% elseif strcmpi(shift_est_method, 'optflow')
%     for t1 = 1:T
%         Y_t1 = Y(:,:,t1);
%         D_Y(t1,t1,:) = 0;
%         for t2 = [max(1,t1-max_shift):t1-1 t1+1:min(T,t1+max_shift)]
%             Y_t2 = Y(:,:,t2);
%             D_Y(t1, t2, :) = optflow_2D(Y_t1, Y_t2, c_ofsigma*abs(t2-t1));
%         end
%     end
% end
% toc;
% Compute errors in the estimates.
% ED_Y = sqrt(sum((D_Y - D_true) .^ 2, 3));

%% Generate SPAD samples
conlog('Generating SPAD samples');
R = 1;
B = false([L T R]);
for t = 1:T
    Yt = Y(:,:,t);
    for r = 1:R
        B(:,:,t,r) = get_spad_shot(Yt);
    end
end

% Ground truth reconstructions.
Yhat_noalign = radiance_mle_spad(align_and_average_2D(B(:,:,:,1), zeros(T, 2), @mean), T, T);
Yhat_gt = radiance_mle_spad(align_and_average_2D(B(:,:,:,1), D0, @mean), T, T);

%% Set up different smoothing windows
conlog('Estimating frame displacements at multiple time scales');
tic;
max_sigma = max(L);
c_ofsigma = 0.25;
enforce_vsc = true;
enforce_lowrank = false;
BS = [1];
Y_mle = zeros([L T R numel(BS)]);
D_Ymle = nan(T, T, 2, R, numel(BS));
D_Ymle_vsc = nan(size(D_Ymle));
for i_bs = 1:numel(BS)
    mle_bsize = BS(i_bs);
    conlog(sprintf('mle_bsize %d', mle_bsize));
    if mle_bsize > 1
        B_bavg = moving_average_1D(double(B), 3, mle_bsize, 1, 'mean');
        Y_mle(:,:,mle_bsize:end,:,i_bs) = radiance_mle_spad(B_bavg, mle_bsize, mle_bsize);
    else
        Y_mle(:,:,:,:,i_bs) = double(B);
    end
    
    if strcmpi(shift_est_method, 'phasecorr')
        Y_mle_F = zeros(size(Y_mle));
        for t = mle_bsize:T
            for r = 1:R
                Y_mle_F(:,:,t,r,i_bs) = dft_2D(Y_mle(:,:,t,r,i_bs), dft_window_fn);
            end
        end
    end
    
    for t1 = mle_bsize:T
        D_Ymle(t1,t1,:,i_bs) = 0;
        for t2 = [max(mle_bsize,t1-half_shift_window):t1-1 t1+1:min(T,t1+half_shift_window)]
            if strcmpi(shift_est_method, 'phasecorr')
                Y_mle_F_t1 = Y_mle_F(:,:,t1,:,i_bs);
                Y_mle_F_t2 = Y_mle_F(:,:,t2,:,i_bs);
                for r = 1:R
                    est_d = phasecorr_dft_2D(Y_mle_F_t1(:,:,1,r), ...
                                             Y_mle_F_t2(:,:,1,r), ...
                                             'plane_fit');
                    D_Ymle(t1, t2, :, r, i_bs) = reshape(est_d, [1 1 2]);
                end
            else
                Y_mle_t1 = Y_mle(:,:,t1,:,i_bs);
                Y_mle_t2 = Y_mle(:,:,t2,:,i_bs);
                for r = 1:R
                    if strcmpi(shift_est_method, 'optflow')
                        est_d = optflow_2D(Y_mle_t1(:,:,1,r), ...
                                           Y_mle_t2(:,:,1,r), ...
                                           min(max_sigma, c_ofsigma*abs(t2-t1)));
                        D_Ymle(t1, t2, :, r, i_bs) = reshape(est_d, [1 1 2]);
                    end
                end
            end
        end
    end
    
    if enforce_vsc
        for r = 1:R
            D_Ymle_vsc(:,:,:,r,i_bs) = force_vectorsum_equality(D_Ymle(:,:,:,r,i_bs), 1);
        end
    end
    
    if enforce_lowrank
        for r = 1:R
            D_Ymle(mle_bsize+1:end, mle_bsize+1:end, 1, r, i_bs) = ...
                        lowrank(D_Ymle(mle_bsize+1:end, mle_bsize+1:end, 1, r, i_bs), ...
                                rank(D_true(:,:,1)));
            D_Ymle(mle_bsize+1:end, mle_bsize+1:end, 2, r, i_bs) = ...
                        lowrank(D_Ymle(mle_bsize+1:end, mle_bsize+1:end, 2, r, i_bs), ...
                                rank(D_true(:,:,2)));
        end
    end
end
toc;

% Compute errors in the estimates.
ED_Ymle = sqrt(sum((D_Ymle - D_true(:,:,:,ones(R,1),ones(numel(BS),1))) .^ 2, 3));
ED_Ymle_vsc = sqrt(sum((D_Ymle_vsc - D_true(:,:,:,ones(R,1),ones(numel(BS),1))) .^ 2, 3));

%% Combine the observations into a single image
conlog('Aligning observed images according to estimated shifts');
tic;
Dhat = zeros(T, 2, R, numel(BS));
Dhat_vsc = zeros(T, 2, R, numel(BS));
for i_bs = 1:numel(BS)
    for r = 1:R
        Dhat(:,:,r,i_bs) = estimated_locs_2D(D_Ymle(:,:,:,r,i_bs));
        Dhat_vsc(:,:,r,i_bs) = estimated_locs_2D(D_Ymle_vsc(:,:,:,r,i_bs));
    end
end

Yhat = zeros([L numel(BS)]);
Yhat_vsc = zeros([L numel(BS)]);
for i_bs = 1:numel(BS)
    valid_idxs = ~any(isnan(Dhat(:,:,1,i_bs)), 2);
    Yhat(:,:,i_bs) = radiance_mle_spad(align_and_average_2D(B(:,:,valid_idxs,1), Dhat(valid_idxs,:,1,i_bs)), ...
                                       T, T);
    if enforce_vsc
        Yhat_vsc(:,:,i_bs) = radiance_mle_spad(align_and_average_2D(B(:,:,valid_idxs,1), Dhat_vsc(valid_idxs,:,1,i_bs)), ...
                                               T, T);
    end
end
toc;

conlog('Done!');