clearvars;
close all;

dft_window_fn = @box_window_1D;
shift_est_method = 'phasecorr';
L = 256;
input_fn = @sin_1D;
offset = 0;
if isequal(input_fn, @step_1D)
    offset = L/2;
end

brightness = 0.2;
contrast = 0.8;

vs = 0.01 * 5; % velocity of the scene
T = 512;

%% Generate the true shifts
D0 = offset - vs + cumsum(vs(ones(T, 1)));
[t2, t1] = meshgrid(1:T);
D_true = D0(t2) - D0(t1);

%% Generate input signal.
Y = zeros(L, T);
for t = 1:T
    Y(:,t) = input_fn(L, D0(t), brightness, contrast);
end

%% Generate baseline shift matrix
max_shift = vs * T;
c_ofsigma = vs;
D_Y = zeros(T, T);
for t1 = 1:T
    for t2 = 1:T
        if strcmpi(shift_est_method, 'optflow')
            D_Y(t1, t2) = optflow_1D(Y(:,t1), Y(:,t2), min(max_shift, c_ofsigma*abs(t2-t1)), 'wmean');
        elseif strcmpi(shift_est_method, 'phasecorr')
            D_Y(t1, t2) = phasecorr_1D(Y(:,t1), Y(:,t2), 'plane_fit', dft_window_fn);
        end
    end
end

%% Generate SPAD samples from input signal.
R = 1;
B = false(L, T, R);
for t = 1:T
    Yt = Y(:,t);
    for r = 1:R
        B(:,t,r) = get_spad_shot(Yt);
    end
end

%% Set up different smoothing windows
enforce_vsc = true;
enforce_lowrank = false;
BS = [1 2 4 8 16 32 64 128];
Y_mle = zeros(L, T, R, numel(BS));
D_Ymle = nan(T, T, R, numel(BS));
for i_bs = 1:numel(BS)
    mle_bsize = BS(i_bs);
    fprintf('%s: mle_bsize %d\n', datestr(datetime('now')), mle_bsize);
    if mle_bsize > 1
        B_bavg = moving_average_1D(double(B), 2, mle_bsize, 1, 'mean');
        Y_mle(:,mle_bsize:end,:,i_bs) = radiance_mle_spad(B_bavg, mle_bsize, mle_bsize);
    else
        Y_mle(:,:,:,i_bs) = double(B);
    end
    
    for t1 = mle_bsize:T
        D_Ymle(t1,t1,:,i_bs) = 0;
        for t2 = [mle_bsize:t1-1 t1+1:T]
            Y_mle_t1 = Y_mle(:,t1,:,i_bs);
            Y_mle_t2 = Y_mle(:,t2,:,i_bs);
            for r = 1:R
                if strcmpi(shift_est_method, 'phasecorr')
                    D_Ymle(t1, t2, r, i_bs) = phasecorr_1D(Y_mle_t1(:,1,r), ...
                                                           Y_mle_t2(:,1,r), ...
                                                           'plane_fit', dft_window_fn);
                elseif strcmpi(shift_est_method, 'optflow')
                    D_Ymle(t1, t2, r, i_bs) = optflow_1D(Y_mle_t1(:,1,r), ...
                                                         Y_mle_t2(:,1,r), ...
                                                         c_ofsigma*abs(t2-t1), 'wmean');
                end
            end
        end
    end
    
    if enforce_vsc
        for r = 1:R
            D_Ymle(:,:,r,i_bs) = force_vectorsum_equality(D_Ymle(:,:,r,i_bs), 1);
        end
    end
    
    if enforce_lowrank
        for r = 1:R
            D_Ymle(mle_bsize+1:end, mle_bsize+1:end, r, i_bs) = ...
                        lowrank(D_Ymle(mle_bsize+1:end, mle_bsize+1:end, r, i_bs), ...
                                rank(D_true));
        end
    end
end

%% Calculate the relative error in the estimates
AE_DY = (D_Y - D_true);
AE_DYmle = (D_Ymle - D_true(:,:,ones(R,1),ones(numel(BS),1)));
RE_DY = AE_DY ./ abs(1e-10 + D_true);
RE_DYmle = AE_DYmle ./ abs(1e-10 + D_true(:,:,ones(R,1),ones(numel(BS),1)));

%% Show the true and the estimated shifts
fig_height = 3*T;
fig_width = 2*T;
figure('Position', [0 0 fig_width fig_height]);
colormap(gray);
sp1 = subplot(3, 4, 2);
image(Y', 'CDataMapping', 'scaled');
xlabel('Spatial index');
ylabel('Time');
title('Radiance');
colorbar;
sp2 = subplot(3, 4, [3 4]);
image(D_true, 'CDataMapping', 'scaled');
ylabel('Source frame');
xlabel('Target frame');
title(sprintf('True shift matrix'));
colorbar;
sp3 = subplot(3, 4, 6);
image(Y', 'CDataMapping', 'scaled');
xlabel('Spatial index');
ylabel('Time');
title('Radiance');
colorbar;
sp4 = subplot(3, 4, [7 8]);
image(D_Y, 'CDataMapping', 'scaled');
ylabel('Source frame');
xlabel('Target frame');
title(sprintf('Baseline shift matrix'));
colorbar;
sp5 = subplot(3, 4, [9 10]);
image(AE_DY, 'CDataMapping', 'scaled');
ylabel('Source frame');
xlabel('Target frame');
title(sprintf('Original error in shift estimates'));
colorbar;
sp6 = subplot(3, 4, [11 12]);
image(RE_DY, 'CDataMapping', 'scaled');
ylabel('Source frame');
xlabel('Target frame');
title(sprintf('Relative error in shift estimates'));
colorbar;
linkaxes([sp1 sp2 sp3 sp4 sp5 sp6], 'y');
drawnow;
print('D_Y', '-dpng');
close();
for i_bs = 1:numel(BS)
    mle_bsize = BS(i_bs);
    figure('Position', [0 0 fig_width fig_height]);
    colormap(gray);
    sp1 = subplot(3, 4, 2);
    image(Y', 'CDataMapping', 'scaled');
    xlabel('Spatial index');
    ylabel('Time');
    title('Radiance');
    colorbar;
    sp2 = subplot(3, 4, [3 4]);
    image(D_true, 'CDataMapping', 'scaled');
    ylabel('Source frame');
    xlabel('Target frame');
    title(sprintf('True shift matrix'));
    colorbar;
    sp3 = subplot(3, 4, 6);
    image(Y_mle(:,:,1,i_bs)', 'CDataMapping', 'scaled');
    xlabel('Spatial index');
    ylabel('Time');
    title('MLEs');
    colorbar;
    sp4 = subplot(3, 4, [7 8]);
    image(D_Ymle(:,:,1,i_bs), 'CDataMapping', 'scaled');
    ylabel('Source frame');
    xlabel('Target frame');
    title(sprintf('Recovered shift matrix for MLE block size %d', mle_bsize));
    colorbar;
    sp5 = subplot(3, 4, [9 10]);
    image(AE_DYmle(:,:,1,i_bs), 'CDataMapping', 'scaled');
    ylabel('Source frame');
    xlabel('Target frame');
    title(sprintf('Original error in shift estimates'));
    colorbar;
    sp6 = subplot(3, 4, [11 12]);
    image(RE_DYmle(:,:,1,i_bs), 'CDataMapping', 'scaled');
    ylabel('Source frame');
    xlabel('Target frame');
    title(sprintf('Relative error in shift estimates'));
    colorbar;
    linkaxes([sp1 sp2 sp3 sp4 sp5 sp6], 'y');
    drawnow;
    print(sprintf('D_Ymle_B%d', mle_bsize), '-dpng');
    close();
end

%% End
conlog('Done!');