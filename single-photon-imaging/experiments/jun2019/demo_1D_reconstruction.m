clearvars; close all;

input_types = {'sin', 'step', 'smooth', 'power_law_fft'};
mle_block_sizes   = [1 2 4 8 16 32 64];
mle_block_strides = [1 1 1 1  1  1  1];
assert(numel(mle_block_sizes) == numel(mle_block_strides));
M = numel(input_types);
N = numel(mle_block_sizes);

bias_ymle = zeros(N, M);
bias_b_from_ymle = zeros(N, M);
iqr_ymle = zeros(N, M);
iqr_b_from_ymle = zeros(N, M);

for i = 1:N
    for j = 1:M
        disp({mle_block_sizes(i) mle_block_strides(i) input_types{j}});
        [by, iy, bb, ib] = simulate_1D_reconstruction(input_types{j}, 256, 256, 64,...
                                                      0.01, 0.2, 0.8, 1, ...
                                                      32, 16, @raised_cosine_window, 'plane_fit',...
                                                      mle_block_sizes(i), mle_block_strides(i));
        bias_ymle(i, j) = by;
        iqr_ymle(i, j)  = iy;
        bias_b_from_ymle(i, j) = bb;
        iqr_b_from_ymle(i, j)  = ib;
    end                                                  
end    