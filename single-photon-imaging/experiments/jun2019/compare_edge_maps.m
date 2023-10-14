clear all; clc;
plots_dir = 'plots';

N_takes = 12;
Ns = int32(logspace(0, (N_takes-1) * log10(2), N_takes));

exposure_time = 1.0;

img_paths = {'greens.jpg', 'cameraman.tif', 'peppers.png', 'football.jpg',...
             'bag.png', 'mri.tif', 'circles.png', 'coins.png', 'canoe.tif'};
for k = 1:numel(img_paths)
    img_path = img_paths{k};
    disp(img_path);
    
    X = im2double(imread(img_path));
    if size(X, 3) > 1
        X = rgb2gray(X);
    end
    true_edge_map = edge(X, 'Canny');

    rmses = zeros(N_takes, 1);
    edge_ious = zeros(N_takes,1);
    X_mle = zeros(size(X, 1), size(X, 2), N_takes);
    edge_maps = false(size(X, 1), size(X, 2), N_takes);
    for i = 1:N_takes
        n = Ns(i)
        snaps = get_binary_snapshots(X, n, exposure_time);
        X_mle(:,:,i) = reconstruct_flux_from_binary_snapshots_mle(snaps, exposure_time);
        d = X_mle(:,:,i) - X;
        rmses(i) = std(d(:));
        edge_maps(:,:,i) = edge(X_mle(:,:,i), 'Canny');
        edge_ious(i) = get_binary_img_IOU(true_edge_map, edge_maps(:,:,i));
    end

    fig = figure;
    semilogx(Ns, rmses, Ns, edge_ious);
    legend('RMSE', 'Edge IOU');
    xlabel('Number of binary snapshots');
    ylabel('RMSE/Edge IOU');
    title(sprintf('MLE Performance for %s', img_path));
    savefig(fig, sprintf(fullfile(plots_dir, sprintf('mle_%s.fig', img_path))));
end