clearvars;
close all;

%% Set up the experiment parameters
dft_window_fn = @box_window_1D;

L = 256;
brightness = 0.2;
contrast = 0.8;
input_fn = @sin_1D;

T = 256;
vs = 0.05;
initial_offset = 0;
if isequal(input_fn, @step_1D)
    initial_offset = L/2;
end

% Set up the ground truth shifts
d0 = [initial_offset; cumsum(vs(ones(T-1, 1)))];
[t2, t1] = meshgrid(1:T);
D_true = d0(t2) - d0(t1);

%% Generate the input and estimate shifts directly from it.
conlog('Generating radiance');
Y = zeros(L, T);
YF = zeros(L, T);
for t = 1:T
    Y(:, t) = input_fn(L, d0(t), brightness, contrast);
    YF(:, t) = dft_1D(Y(:,t), dft_window_fn);
end

dist_fn_of = @(ys1, ys2) optflow_1D(ys1, ys2, 0, 'mean');
dist_fn_pc = @(yF1, yF2) phasecorr_dft_1D(yF1, yF2, 'plane_fit');
D_of_Y = distmat_allpairs_1D(Y, dist_fn_of);
D_pc_Y = distmat_allpairs_1D(YF, dist_fn_pc);

err_DofY = D_of_Y - D_true;
err_DpcY = D_pc_Y - D_true;
dh_ofY = cumsum([initial_offset; diag(D_of_Y, 1)]);
dh_pcY = cumsum([initial_offset; diag(D_pc_Y, 1)]);

%% Generate SPAD samples and estimate shifts from them now.
max_shift = L/2;
c_ofsigma = 0.25;
sigma = c_ofsigma * max_shift;
gk = fspecial('gaussian', [ceil(6*sigma+1) 1], sigma);

conlog('Generating SPAD samples');
R = 1;
B = false(L, T, R);
BF = zeros(L, T, R);
Bs = zeros(L, T, R);
for t = 1:T
    for r = 1:R
        B(:,t,r) = get_spad_shot(Y(:,t));
        BF(:,t,r) = dft_1D(B(:,t,r), dft_window_fn);
        Bs(:,t,r) = conv(double(B(:,t,r)), gk, 'same');
    end
end

dist_fn_of = @(bs1, bs2) optflow_1D(bs1, bs2, 0, 'mean');
D_of_B = distmat_allpairs_1D(Bs, dist_fn_of);
D_pc_B = distmat_allpairs_1D(BF, dist_fn_pc);

err_DofB = D_of_B - D_true(:,:,ones(R,1));
err_DpcB = D_pc_B - D_true(:,:,ones(R,1));
dh_ofB = zeros(T, R);
dh_pcB = zeros(T, R);
for r = 1:R
    dh_ofB(:,r) = cumsum([initial_offset; diag(D_of_B(:,:,r), 1)]);
    dh_pcB(:,r) = cumsum([initial_offset; diag(D_pc_B(:,:,r), 1)]);
end

D_of_B_vsc = zeros(size(D_of_B));
D_pc_B_vsc = zeros(size(D_pc_B));
for r = 1:R
    D_of_B_vsc(:,:,r) = force_vectorsum_equality(D_of_B(:,:,r));
    D_pc_B_vsc(:,:,r) = force_vectorsum_equality(D_pc_B(:,:,r));
end
err_DofBvsc = D_of_B_vsc - D_true(:,:,ones(R,1));
err_DpcBvsc = D_pc_B_vsc - D_true(:,:,ones(R,1));
dh_ofBvsc = zeros(T, R);
dh_pcBvsc = zeros(T, R);
for r = 1:R
    dh_ofBvsc(:,r) = cumsum([initial_offset; diag(D_of_B_vsc(:,:,r), 1)]);
    dh_pcBvsc(:,r) = cumsum([initial_offset; diag(D_pc_B_vsc(:,:,r), 1)]);
end

%% Blockwise averaging
conlog('Blockwise averaging');
S = 7;
BS = 2 .^ (1:S);
Ymle = zeros(L, T, R, S);
Ymle_s = zeros(L, T, R, S);
YmleF = zeros(size(Ymle));
bavg_fn = @mean;
for s = 1:S
    mle_bsize = BS(s);
    conlog(sprintf('MLE block size: %d', mle_bsize));
    if mle_bsize > 1
        for t = 1:mle_bsize-1
            Ymle(:,t,:,s) = bavg_fn(double(B(:,1:t,:)), 2);
        end
        B_bavg = moving_average_1D(double(B), 2, mle_bsize, 1, bavg_fn);
        Ymle(:,mle_bsize:end,:,s) = radiance_mle_spad(B_bavg, mle_bsize, mle_bsize);
    else
        Ymle(:,:,:,s) = double(B);
    end

    for t = 1:T
        for r = 1:R
            YmleF(:,t,r,s) = dft_1D(Ymle(:,t,r,s), dft_window_fn);
            Ymle_s(:,t,r,s) = conv(Ymle(:,t,r,s), gk, 'same');
        end
    end
end

dist_fn_of = @(y1, y2) optflow_1D(y1, y2, 0, 'mean');
D_of_Ymle = distmat_allpairs_1D(Ymle_s, dist_fn_of);
D_pc_Ymle = distmat_allpairs_1D(YmleF, dist_fn_pc);

err_DofYmle = D_of_Ymle - D_true(:,:,ones(R,1),ones(S,1));
err_DpcYmle = D_pc_Ymle - D_true(:,:,ones(R,1),ones(S,1));
dh_ofYmle = zeros(T, R, S);
dh_pcYmle = zeros(T, R, S);
for s = 1:S
    for r = 1:R
        dh_ofYmle(:,r,s) = cumsum([initial_offset; diag(D_of_Ymle(:,:,r,s), 1)]);
        dh_pcYmle(:,r,s) = cumsum([initial_offset; diag(D_pc_Ymle(:,:,r,s), 1)]);
    end
end

D_of_Ymle_vsc = zeros(size(D_of_Ymle));
D_pc_Ymle_vsc = zeros(size(D_pc_Ymle));
for s = 1:S
    for r = 1:R
        D_of_Ymle_vsc(:,:,r,s) = force_vectorsum_equality(D_of_Ymle(:,:,r,s));
        D_pc_Ymle_vsc(:,:,r,s) = force_vectorsum_equality(D_pc_Ymle(:,:,r,s));
    end
end
err_DofYmlevsc = D_of_Ymle_vsc - D_true(:,:,ones(R,1),ones(S,1));
err_DpcYmlevsc = D_pc_Ymle_vsc - D_true(:,:,ones(R,1),ones(S,1));
dh_ofYmlevsc = zeros(T, R, S);
dh_pcYmlevsc = zeros(T, R, S);
for s = 1:S
    for r = 1:R
        dh_ofYmlevsc(:,r,s) = cumsum([initial_offset; diag(D_of_Ymle_vsc(:,:,r,s), 1)]);
        dh_pcYmlevsc(:,r,s) = cumsum([initial_offset; diag(D_pc_Ymle_vsc(:,:,r,s), 1)]);
    end
end

conlog('Done!');