clearvars;
close all;

L = [256 256];
brightness = 0.2;
contrast = 0.8;

vs = [0.1 0.05];
T = 256;

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
end

Y_noalign = align_and_average_2D(Y, zeros(T, 2), @mean);
Y_aligned = align_and_average_2D(Y, D0, @mean);

%% Generate SPAD samples
conlog('Generating SPAD samples');
B = false([L T]);
for t = 1:T
    B(:,:,t) = get_spad_shot(Y(:,:,t));
end

% Ground truth reconstructions.
Yhat_noalign = radiance_mle_spad(align_and_average_2D(B, zeros(T, 2), @mean), T, T);
Yhat_noalign_filtered = wiener_denoise_s(Yhat_noalign, 1, 4);
Yhat_gt = radiance_mle_spad(align_and_average_2D(B, D0, @mean), T, T);
Yhat_noalign_tdenoised = radiance_mle_spad(wiener_denoise_s(wiener_denoise_t0(double(B))), T, T);

%% Set up representations at different levels of smoothing.
conlog('Creating a stack of smoothed versions of images');
tic;
sigmas = [0 2 .^ (-2:4)];
S = numel(sigmas);
sp_kernels = {};
for s = 1:S
    if sigmas(s) > 0
        W = fspecial('gaussian', 2*ceil(2*sigmas(s))+1, sigmas(s));
    else
        W = 1;
    end
    sp_kernels = [sp_kernels W];
end

tw_sizes = 2 .^ [0 1 2 5];
M = numel(tw_sizes);

conlog('Time averaging...');
dB = double(B);
B_bavgs = zeros([L T M]);
for m = 1:M
    conlog(sprintf('wsize = %d', tw_sizes(m)));
    if tw_sizes(m) > 1
        B_bavgs(:,:,:,m) = movmean(dB, [tw_sizes(m)-1 0], 3);
    else
        B_bavgs(:,:,:,m) = dB;
    end
end

conlog('Space averaging...');
sigmas = [0 2 .^ (-2:2:4)];
S = numel(sigmas);
B_s = B_bavgs(:,:,:,:,ones(S,1));
for s = 1:S
    conlog(sprintf('s = %.2f', s));
    for m = 1:M
        conlog(sprintf('wsize = %d', tw_sizes(m)));
        if sigmas(s) > 0
            B_s(:,:,:,m,s) = imgaussfilt(B_bavgs(:,:,:,m), sigmas(s));
        else
            B_s(:,:,:,m,s) = B_bavgs(:,:,:,m);
        end
    end
end
toc;

%% Compute image gradients
conlog('Gradients...');
tic;
B_s_flat = reshape(B_s, [L T*M*S]);
B_s_gradx = reshape(imfilter(B_s_flat, [0, -1, 1], 'replicate'), [L T M S]);
B_s_grady = reshape(imfilter(B_s_flat, [0; -1; 1], 'replicate'), [L T M S]);
toc;

%% Solve optical flow equation with these gradients for all pairs (t1, t2)
conlog('Find shifts at each scale in time/space');
tic;
dt_max = 32;
enforce_vsc = false;
D = zeros(T, T, 2, M, S);
D_vsc = nan(size(D));
for m = 1:M
    conlog(sprintf('wsize = %d', tw_sizes(m)));
    for s = 1:S
        conlog(sprintf('sigma = %f', sigmas(s)));
        I_ms = B_s(:,:,:,m,s);
        Ix_ms = B_s_gradx(:,:,:,m,s);
        Iy_ms = B_s_grady(:,:,:,m,s);
        D_ms = nan(T, T, 2);
        for t1 = 1:T
            D_ms(t1,t1,:) = 0;
            for t2 = [max(1, t1-dt_max):t1-1 t1+1:min(T, t1+dt_max)]
                It = I_ms(:,:,t2) - I_ms(:,:,t1);
                Ix = Ix_ms(:,:,t1);
                Iy = Iy_ms(:,:,t1);
                D_ms(t1,t2,:) = reshape(optflow_2D_direct(Ix, Iy, It),...
                                        [1 1 2]);
            end
        end
        D(:,:,:,m,s) = D_ms;
        if enforce_vsc
            D_vsc(:,:,:,m,s) = force_vectorsum_equality(D_ms, 1);
        end
    end
end

err_D = sqrt(sum((D - D_true(:,:,:,ones(M,1),ones(S,1))) .^ 2, 3));
toc;

%% Combine the observations into a single image
conlog('Aligning observed images according to estimated shifts');
tic;
Dhat = zeros(T, 2, M, S);
Dhat_vsc = zeros(T, 2, M, S);
for m = 1:M
    for s = 1:S
        Dhat(:,:,m,s) = estimated_locs_2D(D(:,:,:,m,s));
        if enforce_vsc
            Dhat_vsc(:,:,m,s) = estimated_locs_2D(D_vsc(:,:,:,m,s));
        end
    end
end

Yhat = zeros([L M S]);
Yhat_vsc = zeros([L M S]);
for m = 1:M
    for s = 1:S
        valid_idxs = ~any(isnan(Dhat(:,:,m,s)), 2);
        Yhat(:,:,m,s) = radiance_mle_spad(align_and_average_2D(B(:,:,valid_idxs), Dhat(valid_idxs,:,m,s)), ...
                                           T, T);
        if enforce_vsc
            Yhat_vsc(:,:,m,s) = radiance_mle_spad(align_and_average_2D(B(:,:,valid_idxs), Dhat_vsc(valid_idxs,:,m,s)), ...
                                                  T, T);
        end
    end
end
toc;

conlog('Done!');