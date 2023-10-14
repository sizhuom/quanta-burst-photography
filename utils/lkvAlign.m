function uv = lkvAlign(im0, im1, std0, iters)
%LKVALIGN Align two images (2D translation) using Lucas-Kanade with
%variance

images = cat(3, im0, im1);
uv = zeros(1, 1, 2);

for i = 1:iters
    [It, Ix, Iy] = partial_deriv(images, uv);
    A = [Ix(:) Iy(:)] ./ std0(:);
    b = -It(:) ./ std0(:);
    x = A \ b;
    uv = uv + reshape(x, size(uv));
end

uv = reshape(uv, [1 2]);

end

