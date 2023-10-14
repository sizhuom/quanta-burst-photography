function uv = lkAlign(im0, im1, iters, uv0)
%LKALIGN Align two images (2D translation) using Lucas-Kanade

if nargin < 4 || isempty(uv0)
    uv = zeros(1, 1, 2, class(uv0));
else
    uv = reshape(uv0, [1 1 2]);
end

for i = 1:iters
    [It, Ix, Iy] = partial_deriv_patch(im0, im1, uv);
    A = [Ix(:) Iy(:)];
    b = -It(:);
    if rank(A) < 2
        break
    end
    x = A \ b;
    if norm(x) > 1
        x = x / norm(x);
    end
    uv = uv + reshape(x, size(uv));
end

uv = reshape(uv, [1 2]);

end

