clearvars;
close all;

%% Set up the experiment parameters
dft_window_fn = @box_window_1D;

L = 256;
brightness = 0.2;
contrast = 0.8;
input_fn = @plmag_1D;

T = 1;
vs = 0.01;
initial_offset = 0;
if isequal(input_fn, @step_1D)
    initial_offset = L/2;
end

% Set up the ground truth shifts
d0 = [initial_offset; cumsum(vs(ones(T-1, 1)))];
[t2, t1] = meshgrid(1:T);
D_true = d0(t2) - d0(t1);

%% Generate the input.
Y = zeros(L, T);
YF = zeros(size(Y));
for t = 1:T
    Y(:, t) = input_fn(L, d0(t), brightness, contrast);
    YF(:, t) = dft_1D(Y(:,t) - mean(Y(:,t)), dft_window_fn);
end
YF_m = abs(YF);
YF_p = angle(YF);

%% Generate SPAD samples.
R = 256;
B = false(L, T, R);
BF = zeros(L, T, R);
for t = 1:T
    for r = 1:R
        B(:,t,r) = get_spad_shot(Y(:,t));
        BF(:,t,r) = dft_1D(B(:,t,r), dft_window_fn);
    end
end

BF_m = abs(BF);
BF_p = angle(BF);

err_m = BF_m - YF_m(:,:,ones(R,1));
err_p = BF_p - YF_p(:,:,ones(R,1));
fs = fftshift((0:L-1)' - L/2) * (1 / L) * 2 * pi;
figure;
plot(fs, mean(err_m, [2 3]), 'o', fs, YF_m, '--k');
figure;
plot(fs, mean(err_p, [2 3]), 'o', fs, YF_m, '--k');
