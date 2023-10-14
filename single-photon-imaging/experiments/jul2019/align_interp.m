clearvars;
close all;

L = 16;
t = 1:L;

x0 = sin_1D(L, 0, 0.1, 0.9, 1);
b0 = get_spad_shot(x0);

shifts = [0.01 0.1 1 L/4 L/2 L];
bcond = 'constant';
extrap = 'edge_vals';
assert(numel(shifts) == 6);

figure;
for i = 1:6
    shifted_x = translate_1D(x0, shifts(i), bcond, extrap);
    subplot(2, 3, i);
    plot(t, x0, '--', t, shifted_x, 'LineWidth', 2);
    legend('original', 'shifted');
    title(sprintf('shift = %.2f', shifts(i)));
end
suptitle(sprintf('%s shift (rightward)', bcond));

figure;
for i = 1:6
    shifted_b = translate_1D(double(b0), shifts(i), bcond, extrap);
    subplot(2, 3, i);
    plot(t, b0, '--', t, shifted_b, 'LineWidth', 2);
    legend('original', 'shifted');
    title(sprintf('shift = %.2f', shifts(i)));
end
suptitle(sprintf('%s shift (rightward)', bcond));